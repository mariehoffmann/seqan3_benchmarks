/*
 * Version with list of relative anchor gaps (classical anchor gap approach) and
 * block-wise accumulated gap lengths. Operator[] is not constant, but jumps
 * forward block-wise to allocate the block where the virtual address points at.
 * Same for insertion and deletion.
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

#define LOG_LEVEL_AB 0

namespace seqan3 {

    template <typename inner_type, unsigned short block_size=32>
    requires alphabet_concept<ranges::v3::value_type_t<inner_type>>
    //in std namespace and renamed! && random_access_range_concept<inner_type> && sized_range_concept<inner_type>
    struct anchor_blocks
    {

    private:
        using gap_decorator_t   = anchor_blocks;

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
        using location_type     = std::tuple<size_type, size_type, size_type>;

        constexpr anchor_blocks()
        {
            data = std::shared_ptr<data_t>(new data_t{});
        };

        constexpr anchor_blocks(anchor_blocks const &) = default;

        constexpr anchor_blocks & operator=(anchor_blocks const &) = default;

        constexpr anchor_blocks (anchor_blocks && rhs) = default;

        constexpr anchor_blocks & operator=(anchor_blocks && rhs) = default;

        ~anchor_blocks() {data->gap_block_list.clear(); data->gap_sums.clear();};

        constexpr anchor_blocks(inner_type * sequence): data{new data_t{sequence}}
        {
            size_type num_blocks = std::max<size_type>(1, sequence->size()/block_size);
            data->gap_sums = gap_sums_type(num_blocks, 0);
            data->gap_block_list = gap_block_list_t(num_blocks, gap_block_type(0));
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
            if (data->sequence != rhs.data->sequence || this->size() != rhs.size())
                return false;
            for (std::uint64_t i = 0; i < data->gap_block_list.size(); ++i)
            {
                if (data->gap_block_list[i].size() != rhs->data->gap_block_list[i].size())
                    return false;
                std::string g_lhs = gap_block_list_hash(data->gap_block_list[i]);
                std::string g_rhs = gap_block_list_hash(rhs->data->gap_block_list[i]);
                if (g_lhs != g_rhs)
                    return false;
            }
            return true;
        }

        bool operator!=(gap_decorator_t & rhs)   // DONE
        {
            return !(*this == rhs);
        }

        void swap(gap_decorator_t & rhs)         // DONE
        {
            data.swap(rhs.data);
        }

        size_type size() const noexcept             // ok
        {
            if (!data->gap_sums.size()){
                if (!data->sequence)
                    return 0;
                return data->sequence->size();
            }
            size_type gap_acc = std::accumulate(data->gap_sums.begin(), data->gap_sums.end(), 0, std::plus<size_type>());
            //std::cout << "current size = " << data->sequence->size() + gap_acc << std::endl;
            return data->sequence->size() + gap_acc;
        }

        size_type max_size() const                  // DONE
        {
            return inner_type{}.max_size() + data->gap_sums.max_size();
        }

        bool empty() const                          // DONE
        {
            return ((!data->sequence) ? true : data->sequence->empty()) && data->gap_sums.size() == 0;
        }

        // position is virtual, gaps block-wise accumulated
        bool insert_gap(size_type const pos, size_type const size=1)
        {
            assert(pos <= this->size());
            std::cout << "enter insert_gap with (pos, size) = (" << pos << ", " << size << ")" << std::endl;
            // locate lower bounding gap
            location_type location{0, 0, 0};
            locate_gap(pos, location);
            size_type block_id = std::get<0>(location);
            size_type gap_id = std::get<1>(location);
            size_type gap_acc = std::get<2>(location);

            gap_block_type block = data->gap_block_list[std::get<0>(location)];
            // case 1: insert new gap
            // case 1 a) is empty or preceeding gap does not overlap with insertion pos
            if (!block.size() || gap_id == block.size() && block[gap_id-1].first + block[gap_id-1].second + gap_acc < pos)
                block.push_back(gap_t{pos, size});
            // case 1 b) not empty, lower_bound points to upstream gap
            // and check for preceeding gap that may overlap with insertion position
            else if (gap_id != block.size() && block[gap_id].first > pos &&
                (!gap_id || block[gap_id-1].first + block[gap_id-1].second + gap_acc < pos))
            {
                data->gap_block_list[block_id].insert(data->gap_block_list[block_id].begin() + gap_id, gap_t{pos,size});
            }
            // case 2: extend existing gap
            // case 2 a) lower_bound coincides with pos
            else if (block[gap_id].first + gap_acc == pos)
            {
                data->gap_block_list[block_id][gap_id].second += size;
            }
            // case 2 b) preceeding gap gets extended
            else if (!gap_id && block[gap_id-1].first + block[gap_id-1].second + gap_acc >= pos)
            {
                data->gap_block_list[block_id][gap_id-1].second += size;
            }

            else std::cout << "ERROR: should not reach this\n";

            // c) update gap_sum for this block
            data->gap_sums[block_id] += size;
            return true;
        }

        // erase gaps from range aligned_seq[pos1;pos2-1]
        bool erase_gap(size_type const pos1, size_type const pos2)
        {
            assert(pos2 <= this->size());
            // a) locate gap
            location_type location{0, 0, 0};
            locate_gap(pos1, location);
            size_type block_id = std::get<0>(location);
            size_type gap_id = std::get<1>(location);
            size_type gap_acc = std::get<2>(location);

            // there is no gap
            assert(data->gap_block_list[block_id].size() == gap_id);

            // virtual start and end positions of lower bounding gap
            size_type gap_vpos_start = data->gap_block_list[block_id][gap_id].first + gap_acc;
            size_type gap_vpos_end = gap_vpos_start + data->gap_block_list[block_id][gap_id].second;

            // TODO: test erase_gap
            // case:  erase complete gap
            if (pos1 == gap_vpos_start && pos2 == gap_vpos_end)
            {
                std::cout << "case: complete gap erasure in block " << block_id << std::endl;
                data->gap_block_list[block_id].erase(data->gap_block_list[block_id].begin() + gap_id);
            }
            // case: erase gap partially by shortening its length
            else if (pos1 < gap_vpos_start)  // virtual gap start
            {
                std::cout << "case: partial erase of gap " << gap_id << " in block " << block_id << std::endl;
                data->gap_block_list[block_id][gap_id] -= pos2 - pos1;
            }
            else
                std::cout << "ERROR: should not reach this part\n";

            // c) update gap_sum for this block
            data->gap_sums[block_id] -= pos2 - pos1;
            return true;
        }

        constexpr reference operator[](size_type const idx) // const noexcept(noexcept((*host)[pos+n]))
        {
            std::cout << "operator[] with idx = " << idx << std::endl;
            assert(idx < this->size());
            // identify gap block
            // TODO: wasn't there a new Pyton-like feature allowing immediate unpacking?
            location_type location{0, 0, 0};
            locate_gap(idx, location);
            size_type block_id = std::get<0>(location);
            size_type gap_id = std::get<1>(location);
            size_type gap_acc = std::get<2>(location);

            // case 1: no gaps in this block
            if (!data->gap_block_list[block_id].size())
                return (*data->sequence)[idx - gap_acc];
            // case 2: lower_bounding idx coincides with gap start
            else if (data->gap_block_list[block_id][gap_id].first + gap_acc == idx)
                return gap::GAP;
            // case 3: preceeding gap overlaps with idx
            else if (gap_id > 0 && data->gap_block_list[block_id][gap_id-1].first + data->gap_block_list[block_id][gap_id-1].second + gap_acc >= idx)
                return gap::GAP;
            else if (gap_id < data->gap_block_list[block_id][gap_id].first + gap_acc)
                return (*data->sequence)[idx - gap_acc];
            else
                std::cout << "ERROR: should reach this case\n";
            return gap::GAP;
        }

        // helper for benchmark only
        void resize(size_t new_size)
        {
            while (this->size() > new_size)
            {
                if ((value_type)(*this)[new_size-1] == gap::GAP)
                    erase_gap(new_size, new_size+1);
                else
                    this->data->sequence.resize(this->data->sequence.size()-1);
            }
        }

    private:

        // Comparator (lambda, but with restriction of signature)
        bool gap_cmp(const gap_t& lhs, const gap_t& rhs)
        {
            return lhs.first < rhs.first;
        }

        // given a virtual position compute lower bounding gap and return block and gap indices.
        // lower bound is the first element in the gap block that is not less than
        // (i.e. greater or equal to) pos, or end of block if no such element is found
        // return (block_id, gap_id, gap_acc)
        void locate_gap(size_type const pos, location_type &location)
        {
            size_type block_id = 0, gap_id = 0, gap_acc = 0;
            // 1. locate block and accumulate gaps
            auto it = data->gap_sums.begin();
            while (it != data->gap_sums.end() && (gap_acc + (*it) + block_id*block_size) < pos)
            {
                gap_acc += *it;
                ++block_id;
            }
//            gap_acc = std::sum(data->gap_sums.begin(), data->gap_sums.begin() + block_id - 1, 0);

            // accumulate gaps before lower bound within designated block
            auto it2 = std::lower_bound(data->gap_block_list[block_id].begin(), data->gap_block_list[block_id].end(), gap_t{pos-gap_acc, 0},
                [] (const auto& lhs, const auto& rhs) {
                    return lhs.first < rhs.first;
                });
            gap_acc = std::accumulate(data->gap_block_list[block_id].begin(), it2,
                gap_acc, [](size_type acc, gap_t gap){ return acc += gap.second;});
            gap_id = static_cast<size_type>(it2 - data->gap_block_list[block_id].begin());
            //return std::make_tuple<size_type, size_type, size_type>(block_id, gap_id, gap_acc);
            std::get<0>(location) = block_id;
            std::get<1>(location) = gap_id;
            std::get<2>(location) = gap_acc;
        }

        struct data_t
        {
            inner_type * sequence{};
            // block-wise accumulation of gap lengths, size doesn't change after sequence assignment
            std::vector<size_t> gap_sums{};
            // nested vector to store gaps by block, i.e. gap_block_list[block id] = gap_block_list
            gap_block_list_t gap_block_list{};

        };
        std::shared_ptr<data_t> data;
    };

} // namespace seqan3
