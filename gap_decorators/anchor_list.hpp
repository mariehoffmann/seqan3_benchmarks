/*
 Version with list of anchor gaps, i.e. a vector with pairs of virtual start position and gap length.
 Gap lengths are not accumulated, but anchor positions not relative (in contrast to classical anchor gap approach?), but virtual, i.e. on every insertion in the middle, tailing gaps need to be corrected by the new positional offset.

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

#define LOG_LEVEL_AL 0

namespace seqan3 {

    template <typename inner_type>
    requires alphabet_concept<ranges::v3::value_type_t<inner_type>>
    //in std namespace and renamed! && random_access_range_concept<inner_type> && sized_range_concept<inner_type>
    struct anchor_list
    {

    private:
        using gap_decorator_t   = anchor_list;

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
        using gap_list_t        = std::vector<gap_t>;

        constexpr anchor_list()
        {
            data = std::shared_ptr<data_t>(new data_t{});
        };

        constexpr anchor_list(anchor_list const &) = default;

        constexpr anchor_list & operator=(anchor_list const &) = default;

        constexpr anchor_list (anchor_list && rhs) = default;

        constexpr anchor_list & operator=(anchor_list && rhs) = default;

        ~anchor_list() = default;

        constexpr anchor_list(inner_type * sequence): data{new data_t{sequence}} {};

        auto begin() noexcept
        {
            return iterator{*this, 0};
        }

        auto end() noexcept
        {
            return iterator{*this, size()};
        }

        bool operator==(gap_decorator_t & rhs) // DONE
        {
            if (data->sequence != rhs.data->sequence || this->size() != rhs.size())
                return false;
            for (std::uint64_t i = 0; i < data->gap_list.size(); ++i)
                if (data->gap_list[i] != rhs.data->gap_list[i])
                    return false;
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
            if (!data->gap_list.size()){
                if (!data->sequence)
                    return 0;
                return data->sequence->size();
            }
            size_type gap_sum = std::accumulate(std::next(data->gap_list.begin()),
                                                data->gap_list.end(), static_cast<size_type>(data->gap_list[0].second),
                                                [](size_type s, gap_t gap){return s + gap.second;});
            return data->sequence->size() + gap_sum;
        }

        size_type max_size() const                  // DONE
        {
            return inner_type{}.max_size() + data->gap_list.max_size();
        }

        bool empty() const                          // DONE
        {
            return ((!data->sequence) ? true : data->sequence->empty()) && data->gap_list.size() == 0;
        }

        // TODO: segfault
        // position is virtual, gaps not accumulated
        bool insert_gap(size_type const pos, size_type const size=1)
        {
            assert(pos <= this->size());
            if (LOG_LEVEL_AL) std::cout << "Enter insert_gap\n";
            auto it = std::lower_bound(data->gap_list.begin(), data->gap_list.end(), gap_t{pos, 0}, [](gap_t lhs, gap_t rhs) -> bool { return lhs.first < rhs.first;});
            if (LOG_LEVEL_AL) std::cout << "lower_bound is end(): " << (it == data->gap_list.end()) << std::endl;
            auto it_succ = it;
            // case: extend gap head
            if (it != data->gap_list.end() && (*it).first == pos)
            {
                if (LOG_LEVEL_AL) std::cout << "case 1: gap head extension\n";
                (*it).second += size;
                ++it_succ;
            }
            // case: merge with preceeding gap
            else if (it != data->gap_list.begin() && ((*(it-1)).first + (*(it-1)).second) >= pos)
            {
                if (LOG_LEVEL_AL) std::cout << "case 2: gap merge\n";
                (*(it-1)).second += size;
            }

            // insert new gap
            else
            {
                if (LOG_LEVEL_AL) std::cout << "case 3: insert new gap\n";
                if (it == data->gap_list.end())
                {
                    if (LOG_LEVEL_AL) std::cout << "\tsubcase: push_back\n";
                    if (LOG_LEVEL_AL) std::cout << "gap_list.size: " << data->gap_list.size() << std::endl;
                    data->gap_list.push_back(gap_t{pos, size});
                    if (LOG_LEVEL_AL) std::cout << "\tpush_back done\n";
                    // updating iterator to successor
                    it_succ = data->gap_list.end();
                }
                else
                {
                    if (LOG_LEVEL_AL) std::cout << "\tsubcase: insert at it-1\n";
                    if (it == data->gap_list.begin())
                    {
                        data->gap_list.insert(it, gap_t{pos, size});
                        it_succ = data->gap_list.begin() + 1;
                    }
                    else
                        data->gap_list.insert(it-1, gap_t{pos, size});
                }
            }
            // update tailing gaps
            for (; it_succ < data->gap_list.end(); ++it_succ)
            {
                if (LOG_LEVEL_AL) std::cout << "update tail" << std::endl;
                (*it_succ).first += size;
            }
            return true;
        }

        bool erase_gap(size_type const pos)
        {
            if ((value_type)(*this)[pos] != gap::GAP)
                return false;
            return erase_gap(pos, pos+1);
        }

        // precondition: as[pos1:pos2-1] == gap_{pos2-pos1}
        bool erase_gap(size_type const pos1, size_type const pos2)
        {
            assert(pos1 < pos2);
            auto it = std::lower_bound(data->gap_list.begin(), data->gap_list.end(), gap_t{pos1, 0},
                                       [](gap_t gap1, gap_t gap2){return gap1.first < gap2.first;});
            auto it_succ = it;
            // case 1: pos1 is not start of gap, i.e. correct iterator position and shorten existing gap
            if ((*it).first != pos1 && it > data->gap_list.begin() && (((*(it-1)).first + (*(it-1)).second) >= pos2))
               (*--it).second -= pos2 - pos1;
            else
            {
                if (LOG_LEVEL_AL) std::cout << "case: head or complete gap erasure\n";
                // case 2: erase complete gap
                if ((*it).first == pos1 && (*it).first + (*it).second == pos2)
                    data->gap_list.erase(it);
                // case 3: erase head, i.e. left shift gaps (anchor position remains, length shortened)
                else // ((*it).first == pos1)
                    (*it).second -= pos2 - pos1;
                ++it_succ;
            }
            // update tailing gap positions
            for (; it_succ < data->gap_list.end(); ++it_succ)
            {
                (*it_succ).first -= pos2-pos1;
            }
            if (LOG_LEVEL_AL)
            {
                std::cout << "updated gap_list is: ";
                for (auto gap : data->gap_list)
                    std::cout << "(" << gap.first << ", " << gap.second << "), ";
                std::cout << std::endl;
            }
            return true;
        }


        inner_type & get_underlying_sequence() const
        {
            return data->sequence;
        }

        //!\brief Set pointer to ungapped sequence and reset gap vector.
        void set_underlying_sequence(inner_type & sequence) const
        {
            data->sequence = sequence;
            data->gap_list.clear();
        }

        constexpr reference operator[](size_type const idx) // const noexcept(noexcept((*host)[pos+n]))
        {
            assert(idx < size());

            if (!data->gap_list.size()) return value_type((*data->sequence)[idx]);
            if (LOG_LEVEL_AL)
            {
                std::cout << "gaps = ";
                for (gap_t gap : data->gap_list)
                    std::cout << "(" << gap.first << ", " << gap.second << "), ";
                std::cout << "sequence.size = " << data->sequence->size() << std::endl;
                std::cout << std::endl;
                std::cout << "[] with idx = " << idx << std::endl;
            }
            // case1: no gaps before position idx
            if (!data->gap_list.size() || idx < data->gap_list[0].first)
            {
                if (LOG_LEVEL_AL) std::cout << "gaps.size == 0 or first gap pos < idx\n";
                return (*data->sequence)[idx];
            }
            if (LOG_LEVEL_AL) std::cout << "acc gaps ...\n";
            // case 2: compute position in gap or between two gaps
            auto it = data->gap_list.begin();
            // the accumulator points to the first postition AFTER the current gap
            size_type acc = (*it).first + (*it).second;
            size_type gap_acc = (*it).second;
            if (LOG_LEVEL_AL) std::cout << "acc init = " << acc << std::endl;
            ++it;
            // skip forward
            while (it != data->gap_list.end() && idx >= acc + (*it).first - ((*(it-1)).first + (*(it-1)).second))
            {
                // add current gap and underlying sequence offset = current_gap_pos - last_gap_end_pos
                acc += (*it).first - ((*(it-1)).first + (*(it-1)).second) + (*it).second;
                gap_acc += (*it).second;
                ++it;
            }
            if (LOG_LEVEL_AL) std::cout << "acc = " << acc << ", gap_acc = " << gap_acc << "\n";
            if (idx < acc || data->sequence->size() == 0)
            {
                if (LOG_LEVEL_AL) std::cout << "return gap\n";
                return gap::GAP;
            }
            else
                return (*data->sequence)[idx - gap_acc];
        }

        bool resize(size_type new_size)
        {
            if (LOG_LEVEL_AL) std::cout << "start resizing ... ";
            assert(new_size <= this->size());
            for (auto pos = this->size() - 1; pos >= new_size; --pos)
            {
                if ((value_type)(*this)[pos] == gap::GAP)
                    erase_gap(pos);
                else
                    data->sequence->resize(data->sequence->size() - 1);
            }
            if (LOG_LEVEL_AL) std::cout << "... final size = " << this->size() << std::endl;
            return true;
        }

    private:
        struct data_t
        {
            inner_type * sequence{};
            gap_list_t gap_list{};

        };
        std::shared_ptr<data_t> data;
    };

} // namespace seqan3
