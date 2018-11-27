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
            size_type num_blocks = std::max<size_type>(1, sequence->size()/block_size + 1);
            if (LOG_LEVEL_AB) std::cout << "num of blocks: " << num_blocks << std::endl;
            data->gap_sums = gap_sums_type(num_blocks, 0);
            if (LOG_LEVEL_AB) {
                std::cout << "init gap_sums with: ["; for (auto gap_sum : data->gap_sums) std::cout << gap_sum << ", ";
                std::cout << "]" << std::endl;
            }
            data->gap_block_list = gap_block_list_t(num_blocks, gap_block_type(0));
            if (LOG_LEVEL_AB) std::cout << "data->gap_block_list.size = " << data->gap_block_list.size() << std::endl;
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
            //std::cout << "current size = data->sequence->size() + gap_acc " << data->sequence->size()<< " + " << gap_acc<< std::endl;
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

        template<class Tuple>
        decltype(auto) sum_components(Tuple const& tuple) {
            auto sum_them = [](auto const&... e)->decltype(auto) {
                return (e+...);
            };
            return std::apply( sum_them, tuple );
        };
        // position is virtual, gaps block-wise accumulated
        bool insert_gap(size_type const pos, size_type const size=1)
        {
            assert(pos <= this->size());
            if (LOG_LEVEL_AB) std::cout << "enter insert_gap with (pos, size) = (" << pos << ", " << size << ")" << std::endl;
            // locate lower bounding gap
            location_type location{0, 0, 0};
            locate_gap(pos, location);
            size_type block_id = std::get<0>(location);
            size_type gap_id = std::get<1>(location);
            size_type gap_acc = std::get<2>(location);
            int i = 0;
            if (LOG_LEVEL_AB) std::cout << i++ << std::endl;
            if (LOG_LEVEL_AB) std::cout << "block_id = " << block_id << ", gap_id = " << gap_id << ", gap_acc = " << gap_acc << std::endl;
            //gap_block_type block = data->gap_block_list[block_id];

            // case: gap outer extension with gap from preceeding block
            /*if (!block_id && !data->gap_block_list[block_id-1].size() && data->gap_block_list[block_id-1].back().first + data->gap_block_list[block_id-1].back().second + gap_acc == pos)
            {
                std::cout << "\tcase: gap outer extension with gap from preceeding block\n";
                --block_id;
                data->gap_block_list[block_id-1][data->gap_block_list[block_id-1].size()-1].second += size;
            }*/
            // case: gap outer extension with preceeding gap from same block
            // returned gap accumulator contains gap length of preceeding gap and does not need to be added to test for
            if (gap_id > 0 && data->gap_block_list[block_id][gap_id-1].first + gap_acc == pos)
            {
                if (LOG_LEVEL_AB) std::cout << "\tcase: gap outer extension with preceeding gap from same block\n";
                data->gap_block_list[block_id][gap_id-1].second += size;
            }
            // case: insert new gap
            //std::cout << i++ << std::endl;
            else if (!data->gap_block_list[block_id].size() ||
                gap_id == data->gap_block_list[block_id].size() ||
                gap_acc + data->gap_block_list[block_id][gap_id].first > pos)
            {
                if (LOG_LEVEL_AB) std::cout << "\tcase: insert new gap\n";
                // a) current block is empty or there is no coinciding or preceeding gap
                if (!data->gap_block_list[block_id].size() || gap_id == data->gap_block_list[block_id].size())
                {
                    if (LOG_LEVEL_AB) std::cout << "\t\tpush_back (" << pos << ", " << size << ")" << std::endl;
                    data->gap_block_list[block_id].push_back(gap_t{pos-gap_acc, size});
                }
                // b) there is a preceeding gap, insert before
                else
                {
                    if (LOG_LEVEL_AB) std::cout << "\t\tinsert (" << pos << ", " << size << ")" << std::endl;
                    data->gap_block_list[block_id].insert(data->gap_block_list[block_id].begin() + gap_id, gap_t{pos, size});
                }

            }
            // case: gap inner extension
            else if ((gap_acc + data->gap_block_list[block_id][gap_id].first) >= pos
                && (gap_acc + data->gap_block_list[block_id][gap_id].first + data->gap_block_list[block_id][gap_id].second) <= pos)
            {
                if (LOG_LEVEL_AB) std::cout << "\tcase: gap inner extension\n";
                data->gap_block_list[block_id][gap_id].second += size;
            }
            else std::cout << "ERROR: should not reach this\n";

            // update gap_sum for this block
            if (LOG_LEVEL_AB) std::cout << "before update of gap_sum[" << block_id << "] = " << data->gap_sums[block_id] << " with " << size << std::endl;
            data->gap_sums[block_id] += size;
            if (LOG_LEVEL_AB) std::cout << "after update of gap_sum[" << block_id << "] = " << data->gap_sums[block_id] << std::endl;

            if (LOG_LEVEL_AB)
            {
                std::cout << "updated gap_sums: [";
                for (auto gap_sum : data->gap_sums) std::cout << gap_sum << ", ";
                std::cout << "]" << std::endl;
            }
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
            if (LOG_LEVEL_AB) {
                std::cout << "entered erase gap with pos1, pos2 = " << pos1 << ", " << pos2 << ", get location = (" << block_id << ", " << gap_id << ", " << gap_acc << ")" << std::endl;
                std::cout << "gap_sums = [";
                for (auto gap_sum : data->gap_sums)
                    std::cout << gap_sum << ", ";
                std::cout << "\n";
            }
            // there is no gap

            assert(data->gap_block_list[block_id].size() != gap_id);
            //-AAAAACGTTGCA----
            //01234567890123456
            // virtual start and end positions of lower bounding gap
            size_type gap_vpos_start = data->gap_block_list[block_id][gap_id].first + gap_acc;
            size_type gap_vpos_end = gap_vpos_start + data->gap_block_list[block_id][gap_id].second;

            // gap addressed by pos1 is shorter than the range to be erased
            assert(gap_vpos_end >= pos2);

            // case:  erase complete gap
            if (pos1 == gap_vpos_start && pos2 == gap_vpos_end)
            {
                if (LOG_LEVEL_AB) std::cout << "case: complete gap erasure in block " << block_id << std::endl;
                data->gap_block_list[block_id].erase(data->gap_block_list[block_id].begin() + gap_id);
            }
            // case: erase gap partially by shortening its length
            else //    if (pos1 < gap_vpos_start)  // virtual gap start
            {
                if (LOG_LEVEL_AB) std::cout << "case: partial erase of gap " << gap_id << " in block " << block_id << std::endl;
                data->gap_block_list[block_id][gap_id].second -= pos2 - pos1;
            }
            //else    std::cout << "ERROR: should not reach this part\n";

            // finally update gap sum statistics for this block
            if (LOG_LEVEL_AB) std::cout << "gap_sum block before update: " << data->gap_sums[block_id] << ", will substract pos2-pos1 = " << (pos2-pos1) << std::endl;
            data->gap_sums[block_id] -= pos2 - pos1;
            if (LOG_LEVEL_AB) std::cout << "gap_sum block after update: " << data->gap_sums[block_id] << std::endl;
            return true;
        }

        // TODO: rework with changed locate logic
        constexpr reference operator[](size_type const idx) // const noexcept(noexcept((*host)[pos+n]))
        {
            //if (LOG_LEVEL_AB) std::cout << "operator[] with idx = " << idx << std::endl;
            assert(idx < this->size());
            // identify gap block
            // TODO: wasn't there a new Pyton-like feature allowing immediate unpacking?
            location_type location{0, 0, 0};
            locate_gap(idx, location);
            size_type block_id = std::get<0>(location);
            size_type gap_id = std::get<1>(location);
            size_type gap_acc = std::get<2>(location);
            //if (LOG_LEVEL_AB) std::cout << "location[" << idx << "] = (" << std::get<0>(location) << ", " << std::get<1>(location) << ", " << std::get<2>(location) << ")\n";
            /*std::cout << "!data->gap_block_list[block_id].size(): " << (!data->gap_block_list[block_id].size()) << std::endl;
            if ((data->gap_block_list[block_id].size())){
                std::cout << "gap_id ==  data->gap_block_list[block_id].size(): " << (gap_id ==  data->gap_block_list[block_id].size()) << std::endl;
                if (!(gap_id ==  data->gap_block_list[block_id].size())){
                    std::cout << "data->gap_block_list[block_id][gap_id].first + gap_acc > idx: " << (data->gap_block_list[block_id][gap_id].first + gap_acc > idx) << std::endl;
                    std::cout << "data->gap_block_list[block_id][gap_id].first = " << data->gap_block_list[block_id][gap_id].first << ", gap_acc = " << gap_acc << ", idx = " << idx << std::endl;
                }
            }*/
            if (!data->gap_block_list[block_id].size() || gap_id ==  data->gap_block_list[block_id].size() ||
                data->gap_block_list[block_id][gap_id].first + gap_acc > idx)
            {
                //if (LOG_LEVEL_AB) std::cout << "case letter, accessed at idx-gap_acc, i.e. " << idx << " - " << gap_acc << " = " << (idx - gap_acc) << ", current sequence.length = " << data->sequence->size() << std::endl;

                return (*data->sequence)[idx - gap_acc];
            }
            else //if (data->gap_block_list[block_id][gap_id].first + gap_acc <= idx &&
                //data->gap_block_list[block_id][gap_id].first + data->gap_block_list[block_id][gap_id].second + gap_acc > idx)
            {
                //if (LOG_LEVEL_AB) std::cout << "case gap\n";
                return gap::GAP;
            }

            //else
            //    std::cout << "ERROR: should not reach this case\n";
            return gap::GAP;
        }

        // helper for benchmark only
        void resize(size_t new_size)
        {
            if (LOG_LEVEL_AB)
            {
                //std::cout << "enter resize with new_size = " << new_size << ", current size = " << this->size() << std::endl;
            }
            while (this->size() > new_size)
            {
                location_type location;
                locate_gap(new_size, location);
                size_type block_id = std::get<0>(location);
                size_type gap_id = std::get<1>(location);
                //size_type gap_acc = std::get<2>(location); // (value_type)(*this)[i]

                value_type v = (value_type)(*this)[new_size];
                //std::cout << "resize: delete [" << new_size << "] = " << v << std::endl;
                if (v == gap::GAP)
                {
                //    std::cout << "resize: erase_gap at " << new_size << std::endl;

                    erase_gap(new_size, new_size+1);
                }
                else
                {
                //    std::cout << "resize: case resize underlying sequence from " << this->data->sequence->size() << " to " <<  this->data->sequence->size() -1<< std::endl;
                    this->data->sequence->resize(this->data->sequence->size()-1);
                    // shift left all subsequent gaps
                    while (block_id < data->gap_block_list.size())
                    {
                        while (gap_id < data->gap_block_list[block_id].size())
                        {

                            data->gap_block_list[block_id][gap_id].first -= 1;
                            // move to preceeding block
                            if (!gap_id && data->gap_block_list[block_id][gap_id].first % block_size == block_size-1)
                            {
                                size_type gap_len = data->gap_block_list[block_id][gap_id].second;
                                data->gap_block_list[block_id-1].push_back(gap_t(data->gap_block_list[block_id][gap_id].first, gap_len));
                                data->gap_block_list[block_id].erase(data->gap_block_list[block_id].begin());
                                data->gap_sums[block_id] -= gap_len;
                                data->gap_sums[block_id+1] += gap_len;

                            }
                            ++gap_id;
                        }
                        gap_id = 0;
                        ++block_id;
                    }
                }
                if (LOG_LEVEL_AB)
                {
                    std::cout << "resized sequence: has " << data->sequence->size() << " letters and " << std::accumulate(data->gap_sums.begin(), data->gap_sums.end(), 0)<< " gaps \n";
                    for (size_type i = 0; i < this->size(); ++i) std::cout << (value_type)(*this)[i];
                    std::cout << std::endl;
                }
            }
            // remove supernumery tailing blocks
            size_type num_blocks = std::max<size_type>(1, new_size/block_size + 1);
            if (LOG_LEVEL_AB) std::cout << "number of blocks: " << num_blocks << std::endl;
            data->gap_block_list.resize(num_blocks);
            data->gap_sums.resize(num_blocks);
        }

    private:

        // Comparator (lambda, but with restriction of signature)
        bool gap_cmp(const gap_t& lhs, const gap_t& rhs)
        {
            return lhs.first < rhs.first;
        }

    public:
        // given a virtual position compute lower bounding gap and return block and gap indices.
        // behaviour is similar to std::lower_bound, except that when pos is inside a
        // a gap, it returns the earlier starting gap anchor location.
        // If pos does not correspond to an inner gap, the behaviour is identical to
        // std::lower bound and location will be the gap that is not less than
        // (i.e. greater or equal to) pos, or end of block if no such element is found.
        // In addition to the block_id and gap_id the accumulated gap lengths are returned
        // for avoiding repeated computation.
        // location = (block_id, gap_id, gap_acc)
        void locate_gap(size_type const pos, location_type &location)
        {
            size_type block_id = 0, gap_id = 0, gap_acc = data->gap_sums[0];
            gap_t last_gap = (!data->gap_block_list[0].size()) ?  gap_t{0, 0} : data->gap_block_list[0].back();

            // 1. locate block and accumulate gaps
            auto it = data->gap_sums.begin();
//            while (it != data->gap_sums.end() && (gap_acc + last_gap.first + last_gap.second - 1) < pos)
            //
            while (it != data->gap_sums.end() && (((gap_acc + last_gap.first + last_gap.second) < pos + 1) && ((block_id+1)*block_size-1 < pos - gap_acc) ))
            {
                gap_acc += *++it;
                ++block_id;
                last_gap = (!data->gap_block_list[block_id].size()) ? gap_t{0, 0} : data->gap_block_list[block_id].back();
            }
            //std::cout << "\tblock_id = " << block_id << " with size = " << data->gap_block_list[block_id].size() << "\n";
            // reset gap_acc in order to identify true preceeding gap
            if (gap_acc) gap_acc -= data->gap_sums[block_id];
            if (LOG_LEVEL_AB)
            {
                //std::cout << "enter lower_bound with pos -gap_acc, i.e. " << pos << " - " << gap_acc << " = " << pos-gap_acc << std::endl;
                //std::cout << "data->gap_block_list[block_id].size() = " << data->gap_block_list[block_id].size() << std::endl;
                //for (size_type i = 0; i < data->gap_block_list[block_id].size(); ++i)
                //    std::cout << "block contains: (" << data->gap_block_list[block_id][i].first << ", " << data->gap_block_list[block_id][i].second << ")\n";
            }
            // accumulate gaps before lower bound within designated block
            // either do linear search from front to last or binary search with lower_bound, but the need to compute block-local gap_sums
            size_type gap_acc_block = 0;
            auto it2 = data->gap_block_list[block_id].begin();
            if (data->gap_block_list[block_id].size())
            {
                while (it2 != data->gap_block_list[block_id].end() && (gap_acc_block + (*it2).first + (*it2).second - 1 < pos-gap_acc))
                {
                    gap_acc_block += (*it2).second;
                    ++it2;
                }
            }
            //if (LOG_LEVEL_AB) std::cout << "it2 - begin() = " << (it2 - data->gap_block_list[block_id].begin()) << std::endl;
            //std::cout << "it2 == end() ? " << (it2 == data->gap_block_list[block_id].end()) << std::endl;
            /*
            auto it2 = std::lower_bound(data->gap_block_list[block_id].begin(), data->gap_block_list[block_id].end(), gap_t{pos-gap_acc, 0},
                [] (const auto& lhs, const auto& rhs) -> bool {
                    std::cout << "compare: lhs.first = " << lhs.first << ", lhs.second = " << lhs.second << " with " << rhs.first  << ", " << (lhs.first + lhs.second - 1) << " < " << rhs.first << std::endl;
                    return lhs.first + lhs.second - 1 < rhs.first;  // element < value
                });
                */

            //if (LOG_LEVEL_AB) std::cout << "\tit2 for gap_list points to end(): " << (it2 - data->gap_block_list[block_id].end() == data->gap_block_list[block_id].size()) << std::endl;
            gap_acc = std::accumulate(data->gap_block_list[block_id].begin(), it2,
                gap_acc, [](size_type acc, gap_t gap){ return acc += gap.second;});
            gap_id = static_cast<size_type>(it2 - data->gap_block_list[block_id].begin());
            //return std::make_tuple<size_type, size_type, size_type>(block_id, gap_id, gap_acc);
            std::get<0>(location) = block_id;
            std::get<1>(location) = gap_id;
            std::get<2>(location) = gap_acc;
        }
private:
        struct data_t
        {
            inner_type * sequence{};
            // block-wise accumulation of gap lengths, size doesn't change after sequence assignment
            std::vector<size_t> gap_sums{};
            // nested vector to store gaps by block, i.e. gap_block_list[block id] = gap_block_list
            // an eventually tailing gap will be stored in the last block for having consistent gap read behaviour in all blocks
            gap_block_list_t gap_block_list{};
        };
        std::shared_ptr<data_t> data;
    };

} // namespace seqan3
