/*
 * Version with list of relative anchor gaps (classical anchor gap approach) and
 * block-wise accumulated gap lengths. In this version meta blocks store the accumulated
 * gap sum of the previous blocks. Operator[] is not constant, but log(b) + k/b
 * by doing binary search to determine the designated block and than iterates over
 * all contained gaps within that block.
 * On insertion/erasure tailing meta blocks need to be updated: theta(b)
 */

#pragma once

#include <algorithm>
#include <initializer_list>
#include <type_traits>
#include <vector>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/range/all.hpp>
#include <seqan3/range/detail/random_access_iterator.hpp>

#include <sdsl/sd_vector.hpp>
#include <sdsl/util.hpp>

#define LOG_LEVEL_AB2 1

namespace seqan3 {

    template <typename inner_type, unsigned short block_size=32>
    requires alphabet_concept<ranges::v3::value_type_t<inner_type>>
    //in std namespace and renamed! && random_access_range_concept<inner_type> && sized_range_concept<inner_type>
    struct anchor_blocks2
    {

    private:
        using gap_decorator_t   = anchor_blocks2;

    public:
        using alphabet_type = typename ranges::v3::value_type_t<inner_type>;
        using value_type = gapped<alphabet_type>;

        using reference = value_type;

        using const_reference = const reference;
        using iterator = detail::random_access_iterator<gap_decorator_t>;

        using const_iterator = iterator;
        using difference_type = typename ranges::v3::difference_type_t<inner_type>;

        using size_type = typename ranges::v3::size_type_t<inner_type>;

        using gap_t             = std::pair<size_type, size_type>;
        using gap_sums_type     = std::vector<size_type>;
        using gap_block_type    = std::vector<gap_t>;
        using gap_block_list_t  = std::vector<gap_block_type>;

        // store location information, i.e. block and gap id and accumulated gap symbols
        struct location_type
        {
            size_type block_id;
            size_type gap_id;
            size_type gap_acc;
            bool is_gap;
            location_type() : block_id(0), gap_id(0), gap_acc(0), is_gap(false){};
            location_type(size_type _block_id, size_type gap_id, size_type gap_acc, bool is_gap) :
                block_id(block_id), gap_id(gap_id), gap_acc(gap_acc), is_gap(is_gap){};
        };

        constexpr anchor_blocks2()
        {
            size_type num_blocks = 1;
            gap_sums = gap_sums_type(num_blocks, 0);
            gap_block_list = gap_block_list_t(num_blocks, gap_block_type(0));
            //block_list_t b(block_list_t(2, block_t(0)));

        };

        constexpr anchor_blocks2(anchor_blocks2 const &) = default;

        constexpr anchor_blocks2 & operator=(anchor_blocks2 const &) = default;

        constexpr anchor_blocks2 (anchor_blocks2 && rhs) = default;

        constexpr anchor_blocks2 & operator=(anchor_blocks2 && rhs) = default;

        ~anchor_blocks2() {gap_block_list.clear(); gap_sums.clear();};

        constexpr anchor_blocks2(inner_type * sequence): sequence{sequence}
        {
            size_type num_blocks = std::max<size_type>(1, sequence->size()/block_size + 1);
            if (LOG_LEVEL_AB2) std::cout << "num of blocks: " << num_blocks << std::endl;
            gap_sums = gap_sums_type(num_blocks, 0);
            if (LOG_LEVEL_AB2) {
                std::cout << "init gap_sums with: ["; for (auto gap_sum : gap_sums) std::cout << gap_sum << ", ";
                std::cout << "]" << std::endl;
            }
            gap_block_list = gap_block_list_t(num_blocks, gap_block_type(0));

            if (LOG_LEVEL_AB2) std::cout << "gap_block_list.size = " << gap_block_list.size() << std::endl;
        };

        auto begin() noexcept
        {
            return iterator{*this, 0};
        }

        auto end() noexcept
        {
            return iterator{*this, size()};
        }

        std::string gap_block_list_hash(gap_block_list_t & gap_block_list)
        {
            return std::accumulate(std::next(gap_block_list.begin()), gap_block_list.end(),
                    std::to_string(gap_block_list[0].first) + "," + std::to_string(gap_block_list[0].second),
                        [](std::string a, gap_t gap){
                            return a + ";" +  std::to_string(gap.first) + "," + std::to_string(gap.second);
                    });
        }

        bool operator==(gap_decorator_t & rhs)
        {
            if (sequence != rhs.sequence || this->size() != rhs.size())
                return false;
            for (std::uint64_t i = 0; i < gap_block_list.size(); ++i)
            {
                if (gap_block_list[i].size() != rhs->gap_block_list[i].size())
                    return false;
                std::string g_lhs = gap_block_list_hash(gap_block_list[i]);
                std::string g_rhs = gap_block_list_hash(rhs->gap_block_list[i]);
                if (g_lhs != g_rhs)
                    return false;
            }
            return true;
        }

        bool operator!=(gap_decorator_t & rhs)
        {
            return !(*this == rhs);
        }

        void swap(gap_decorator_t & rhs)
        {
            sequence.swap(rhs.sequence);
            gap_sums.swap(rhs.gap_sums);
            gap_block_list.swap(rhs.gap_block_list);
        }

        size_type size() const noexcept
        {
            if (!gap_sums.size()){
                if (!sequence)
                    return 0;
                return sequence->size();
            }
            return sequence->size() + gap_sums.back();
        }

        size_type max_size() const
        {
            return inner_type{}.max_size() + gap_sums.max_size();
        }

        bool empty() const
        {
            return ((!sequence) ? true : sequence->empty()) && gap_sums.size() == 0;
        }

        template<class Tuple>
        decltype(auto) sum_components(Tuple const& tuple) {
            auto sum_them = [](auto const&... e)->decltype(auto) {
                return (e+...);
            };
            return std::apply( sum_them, tuple );
        };
        /* position is virtual, gaps block-wise accumulated
        // cases:    i        ii  iii          iv
        //           |        |   |             |
        //         [  (pos1,len1)  (pos2, len2)  ]
        // i) prepend new gap into designated block
        *       <=> gap accumulator of loc is identical to the one of previous block
        *  ii) extend existing gap
        *       <=> pos is in range of virtual start and end+1 (inclusive!) address of located gap
        *  iii) insert between existing gaps or append
        *       <=> loc.gap_id refers to gap ending before position index
        */
        bool insert_gap(size_type const pos, size_type const size=1)
        {
            if (LOG_LEVEL_AB2) std::cout << "Enter insert_gap with pos = " << pos << std::endl;
            assert(pos <= this->size());
            location_type loc{};
            if (LOG_LEVEL_AB2) std::cout << "h1\n";
            // binary search to map pos to its designated block
            locate_gap(pos, loc);
            if (LOG_LEVEL_AB2) std::cout << "h2\n";
            // get gap accumulator of previous block
            //if (LOG_LEVEL_AB2)
            std::cout << "location for pos = " << pos << ": (" << loc.block_id << ", " << loc.gap_id << ", " << loc.gap_acc << ", " << loc.is_gap << ")\n";
            size_type gap_acc_pred = (loc.block_id) ? gap_sums[loc.block_id-1] : 0;
            std::cout << "num blocks in block_list = " << gap_block_list.size() << std::endl;
            // i) prepend
            if (loc.gap_acc == gap_acc_pred)
            {
                //if (LOG_LEVEL_AB2)
                std::cout << "h2 a1\n";
                // what cast is happening that vector.begin() results in segmentation fault, but + 0 not?
                gap_t gap{pos - loc.gap_acc, size};
                std::cout << "h2 a2\n";
                std::cout << "gap_block_list[loc.block_id].size = " << gap_block_list[loc.block_id].size() << std::endl;

                gap_block_list[loc.block_id].insert(gap_block_list[loc.block_id].begin() + 0, gap);
//                gap_block_list[loc.block_id].insert(gap_block_list[loc.block_id].begin(), gap_t{static_cast<size_type>(pos - loc.gap_acc), size});
                if (LOG_LEVEL_AB2) std::cout << "h2 a3\n";
            }
            // ii) extend
            else if (pos <= gap_block_list[loc.block_id][loc.gap_id].first + loc.gap_acc)
            {
                if (LOG_LEVEL_AB2) std::cout << "h2 b1\n";
                gap_block_list[loc.block_id][loc.gap_id].second += size;
                if (LOG_LEVEL_AB2) std::cout << "h2 b2\n";
            }
            // iii) insert in the middle or append new gap
            else
            {
                //if (LOG_LEVEL_AB2)
                std::cout << "h2 c1\n";
                gap_block_list[loc.block_id].insert(gap_block_list[loc.block_id].begin() + loc.gap_id, gap_t{static_cast<size_type>(pos - loc.gap_acc), size});
                //if (LOG_LEVEL_AB2)
                std::cout << "h2 c2\n";
            }
            if (LOG_LEVEL_AB2) std::cout << "h3\n";
            // update block statistics
            for (size_type i = loc.block_id; i < gap_sums.size();)
                gap_sums[i++] += size;

            if (LOG_LEVEL_AB2)
            {
                std::cout << "updated gap_sums: [";
                for (auto gap_sum : gap_sums) std::cout << gap_sum << ", ";
                std::cout << "]" << std::endl;
                std::cout << "updated gap_block_list: [";
                for (auto gap_block : gap_block_list)
                {
                    std::cout << "[";
                    for (auto gap : gap_block)
                        std::cout << "(" << gap.first << ", " << gap.second << "), ";
                    std::cout << "]";
                }
                std::cout << "]" << std::endl;
            }
            return true;
        }

        // erase gaps from range aligned_seq[pos1;pos2-1]
        // invariant: gap start remains unchanged if not deleted completely
        bool erase_gap(size_type const pos1, size_type const pos2)
        {
            assert(pos2 <= this->size());
            // a) locate gap
            location_type location{0, 0, 0, false};
            locate_gap(pos1, location);
            if (LOG_LEVEL_AB2) std::cout << "location for pos1 = " << pos1 << ": (" << location.block_id << ", " << location.gap_id << ", " << location.gap_acc << ", " << location.is_gap << ")\n";

            // location has to point to existing gap
            assert(location.is_gap);

            // gap range to be deleted has to be consecutive
            assert(pos1 >= gap_block_list[location.block_id][location.gap_id].first + location.gap_acc - gap_block_list[location.block_id][location.gap_id].second);
            assert(pos2 <= gap_block_list[location.block_id][location.gap_id].first + location.gap_acc);

            // to be deleted range has to be inside located gap
            size_type vend = gap_block_list[location.block_id][location.gap_id].first + location.gap_acc + 1;
            size_type vbeg = vend - gap_block_list[location.block_id][location.gap_id].second - 1;

            if (LOG_LEVEL_AB2)
            {
                std::cout << "initial gap list: ";
                for (gap_t gap : gap_block_list[location.block_id])
                    std::cout << "(" << gap.first << ", " << gap.second << "), ";

                std::cout << "\nvstart = " << vbeg << ", vend = " << vend << std::endl;
            }
            assert(pos1 >= vbeg && pos2 <= vend);

            // case 1: delete complete gap
            if (pos1 == vbeg && pos2 == vend)
                gap_block_list[location.block_id].erase(gap_block_list[location.block_id].begin() + location.gap_id);
            // case 2: decrease gap length
            else
                gap_block_list[location.block_id][location.gap_id].second -= pos2 - pos1;
            if (LOG_LEVEL_AB2)
            {
                std::cout << "new gap list: ";
                for (gap_t gap : gap_block_list[location.block_id])
                    std::cout << "(" << gap.first << ", " << gap.second << "), ";
                std::cout << std::endl;
            }
            // update block statistics
            for (size_type i = location.block_id; i < gap_sums.size();)
                gap_sums[i++] -= pos2 - pos1;
            if (LOG_LEVEL_AB2)
            {
                std::cout << "new gap sums: ";
                for (auto gap_sum : gap_sums)
                    std::cout << gap_sum << ", ";
                std::cout << std::endl;
            }
            return true;
        }

        // TODO: rework with changed locate logic
        constexpr reference operator[](size_type const idx) // const noexcept(noexcept((*host)[pos+n]))
        {
            assert(idx < this->size());
            // identify gap block
            location_type loc;
            locate_gap(idx, loc);
            if (LOG_LEVEL_AB2) std::cout << "idx = " << idx << ", loc = (" << loc.block_id << ", " << loc.gap_id << ", " << loc.gap_acc << ", " << loc.is_gap << ")\n";
            if (loc.is_gap)
                return gap::GAP;
            return (*sequence)[idx - loc.gap_acc];
        }

    private:

        // Comparator (lambda, but with restriction of signature)
        bool gap_cmp(const gap_t& lhs, const gap_t& rhs)
        {
            return lhs.first < rhs.first;
        }

        // given a virtual position compute lower bounding gap and return block,
        // and gap indices, and the accumulated gap lengths of previous and
        // enclosing gaps. If pos points into gap, the enclosing gap location (via
        // its indices in the gap_block_list) is returned, otherwise the succeeding
        // gap or end position.
        void locate_gap(size_type const pos, location_type &location)
        {
            // computes the virtual sequence size given a block index and gap prefixes
            std::function<size_type (size_type)> upper_bound = [=](size_type i) { return (i+1)*block_size + this->gap_sums[i]; };

            // a) binary search to identify designated block
            size_type mid = gap_sums.size()/2;
            // number of gaps and alphabet symbols until end of this block
            size_type num_symbols = upper_bound(mid);
            // number of gaps and alphabet symbols until end of preceeding block
            size_type num_symbols_pred = (mid) ? upper_bound(mid-1) : 0;
            while (pos < num_symbols_pred || pos >= num_symbols)
            {
                if (pos < num_symbols_pred && mid)
                    mid /= 2;
                else if (mid < gap_sums.size() - 1)
                    mid = (mid + gap_sums.size())/2;  // TODO: check for floor for fractional numbers
                num_symbols = upper_bound(mid);
                num_symbols_pred = (mid) ? upper_bound(mid-1) : 0;
            }
            // b) locate upper bounding gap within block, if there is no upper bound gap_id points to end
            size_type gap_id = 0;
            size_type gap_acc = (mid) ? gap_sums[mid-1] : 0;
            if (LOG_LEVEL_AB2) std::cout << "gap acc = "  << gap_acc << std::endl;
            // accumulate gaps as long 'pos' refers to position before or inside the current gap
            if (gap_id < gap_block_list[mid].size())
                if (LOG_LEVEL_AB2) std::cout << "pos = " << pos << ", pos after next gap: " << gap_acc + gap_block_list[mid][gap_id].second + gap_block_list[mid][gap_id].first << std::endl;
            while (gap_id < gap_block_list[mid].size() && pos >= gap_acc + gap_block_list[mid][gap_id].first)
            {
                if (LOG_LEVEL_AB2) std::cout << "enter loop with gap_id = " << gap_id << " \n";
                gap_acc += gap_block_list[mid][gap_id].second;
                if (gap_id < gap_block_list[mid].size())
                {
                    if (LOG_LEVEL_AB2) std::cout << "pos >= gap_acc + gap_block_list[mid][gap_id].first ? " <<  (pos >= gap_acc + gap_block_list[mid][gap_id].first) << "\n";
                    if (LOG_LEVEL_AB2) std::cout << "gap_id < gap_block_list[mid].size() && pos >= gap_acc + gap_block_list[mid][gap_id].first ? " << (gap_id < gap_block_list[mid].size() && pos >= gap_acc + gap_block_list[mid][gap_id].first) << "\n";
                }
                if (LOG_LEVEL_AB2) std::cout << "gap_id = " << gap_id << ", gap_acc = " << gap_acc << "\n";
                ++gap_id;

            }
            if (LOG_LEVEL_AB2) std::cout << "gap_id after loop: " << gap_id << std::endl;
            location.block_id = mid;

            location.is_gap = false;
            if (!gap_id && gap_id < gap_block_list[location.block_id].size() && pos >= gap_acc + gap_block_list[location.block_id][gap_id].first)
                location.is_gap = true;
            // increment gap_id if pos refers to position inside gap, i.e.
            if (gap_id && pos >= gap_acc + gap_block_list[mid][gap_id-1].first - gap_block_list[mid][gap_id-1].second && pos < gap_acc + gap_block_list[mid][gap_id-1].first)
            {
                --gap_id;
                location.is_gap = true;
            }
            // revert increment from previous for-loop
//            if (gap_id) --gap_id;
            location.gap_id = gap_id;
            location.gap_acc = gap_acc; //gap_sums[mid-1] : 0;

        }
private:
        inner_type * sequence{};
        // block-wise accumulation of gap lengths, size doesn't change after sequence assignment
        std::vector<size_t> gap_sums{};
        // nested vector to store gaps by block, i.e. gap_block_list[block id] = gap_block_list
        // an eventually tailing gap will be stored in the last block for having consistent gap read behaviour in all blocks
        gap_block_list_t gap_block_list{};

    };

} // namespace seqan3
