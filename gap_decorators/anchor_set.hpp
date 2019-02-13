/*
    Uses an unordered set for the virtual gap positions to retrieve log (k) lookup time.
    Drawback: since the virtual gap positions correspond to the accumulative sum of gap lengths, on update or deletion, the
    worst-case runtime is O(k) for updating subsequent gap positions.
    lower_bound probably does not take advantage of ordering when inserting tuple/pair, because order criterion not obvious.
 */

//#pragma once

#ifndef anchor_set_H

#define anchor_set_H

#include <iostream>
#include <set>
#include <stdlib.h>

#include <range/v3/all.hpp>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>

#define LOG_LEVEL_AS 0
#define _ 0

namespace seqan3 {

// TODO: remove gap_compare, default behaviour is lexicographical comparison
template <typename gap_t>
struct gap_compare {
    bool operator() (const gap_t& lhs, const gap_t& rhs) const {
        return lhs.first < rhs.first;
    }
};

template <typename inner_type>
struct anchor_set
{
public:

    using alphabet_type = typename ranges::v3::value_type_t<inner_type>;
    using value_type = gapped<alphabet_type>;
    using reference = value_type;
    using size_type = typename ranges::v3::size_type_t<inner_type>;
    using gap_t = typename std::pair<size_t, size_t>;
    using set_iterator = typename std::set<gap_t, gap_compare<gap_t>>::iterator;
//    using set_iterator = typename std::set<gap_t>::iterator;

    constexpr anchor_set() = default;

    constexpr anchor_set(anchor_set const &) = default;

    constexpr anchor_set & operator=(anchor_set const &) = default;

    constexpr anchor_set (anchor_set && rhs) = default;

    constexpr anchor_set & operator=(anchor_set && rhs) = default;

    ~anchor_set() = default;

    //constexpr gap_vector_bit(inner_type * sequence): data{new data_t{sequence}}

    constexpr anchor_set(inner_type * sequence): sequence{sequence} {};

    constexpr anchor_set & operator=(inner_type const & sequence)
    {
        sequence = sequence;
        anchors.clear();
    };

    size_type size() const noexcept
    {
        if (anchors.size())
            return (*(anchors.rbegin())).second + sequence->size();
        return sequence->size();
    }

    bool insert_gap(size_type const pos, size_type const size=1)
    {
        if (LOG_LEVEL_AS) std::cout << "called insert with pos = " << pos << ", size = " << size << std::endl;
        if (!size)
            return false;
        //size_type const pos = it - begin();
        assert(pos <= this->size());

        set_iterator it_set = anchors.begin();
        // case 1: extend previous/surrounding gap already existing
        if ((pos < this->size()) && (((value_type)(*this)[pos] == gap::GAP) || (pos > 0 && (*this)[pos-1] == gap::GAP)))
        {
            it_set = anchors.lower_bound(gap_t{pos, 0/*Unused*/});
            if (it_set == anchors.end() || (*it_set).first > pos)
                it_set = std::prev(it_set);
            gap_t gap{(*it_set).first, (*it_set).second + size};
            // merge with successor
            auto it_next = it_set;
            ++it_next;
            if ((*it_set) < (*std::prev(anchors.end())) && (*it_next).first <= (*it_set).first + size - 1)
            {
                // extend gap for *it_next, delete *(it_next+1)
                gap.second += (*it_next).second;
                anchors.erase(it_next);
            }
            anchors.erase(it_set);
            anchors.insert(gap);
        }
        // case 2: create new anchor gap
        else
        {
            gap_t gap{pos, size};
            // pre: pos not in anchor set, find preceeding gap to add accumulated gaps
            if (anchors.size())
            {
                auto it_aux = anchors.lower_bound(gap_t{pos, 0/*Unused*/});
                if (it_aux != anchors.begin())
                    gap.second += (*--it_aux).second;
            }
            anchors.insert(gap);
        }
        // post-processing: reverse update of succeeding gaps
        rupdate(pos, size);
        return true;
    }


    bool erase_gap(size_type const pos)
    {
        return erase_gap(pos, pos+1);
    }

    bool erase_gap(size_type const pos1, size_type const pos2)
    {

        if (LOG_LEVEL_AS) std::cout << "call erase_gap at pos1 = " << pos1 << ", pos2 = " << pos2 << std::endl;
        //size_type pos1 = first - begin(), pos2 = last - begin();
        set_iterator it = anchors.lower_bound(gap_t{pos1, 0/*Unused*/});
        size_type gap_len = get_gap_length(it);

        if (it == anchors.end() || (*it).first > pos1)
            it = std::prev(it);
        // check if [it, it+gap_len[ covers [first, last[
        //if (!(((*it).first <= pos1) && (((*it).first + gap_len) >= pos2))) // [[unlikely]]
        //    throw gap_erase_failure("The range to be erased does not corresponds to a consecutive gap.");
        // case 1: complete gap is deleted
        if (((*it).first == pos1) && (gap_len == pos2-pos1))
            anchors.erase(it);
        // case 2: gap to be deleted in tail or larger than 1 (equiv. to shift tail left, i.e. pos remains unchanged)
        else
        {
            gap_t gap{(*it).first, (*it).second - pos2 + pos1};
            anchors.erase(it);
            anchors.insert(gap);
        }
        // post-processing: forward update of succeeding gaps
        if (LOG_LEVEL_AS) std::cout << "call update ...\n";
        update(pos1, pos2-pos1);
        if (LOG_LEVEL_AS) std::cout << "... update done\n";
        return true;
    }

    void set_underlying_sequence(inner_type * sequence) const
    {
        sequence = sequence;
    }

    constexpr reference operator[](size_type const i) const
    {
        assert(i < size());
        // case 1: no gaps
        if (!anchors.size()) return value_type((*sequence)[i]);
        // case 2: gaps
        set_iterator it = anchors.upper_bound(gap_t{i, _});

        if (it == anchors.begin())
            return value_type((*(sequence))[i]); // since no gaps happen before i

        it = std::prev(it); // here you can be sure that (*it).first <= i by defnition of upper_bound. No need to test.

        size_type gap_len{(*it).second};
        if (it != anchors.begin())
            gap_len -= (*(std::prev(it, 1))).second;

        if (i < (*it).first + gap_len)
            return gap::GAP;
        else
            return value_type((*(sequence))[i - (*it).second]);
        /*
        typename std::set<gap_t, gap_compare<gap_t>>::iterator it = anchors.lower_bound(gap_t{idx, _});
        if (anchors.size() && it == anchors.end())
            it = std::prev(it);

        size_type acc = 0, gap_len = 0;
        if ((*it).first <= idx || (it != anchors.begin() && (*(std::prev(it))).first <= idx))
        {
            if ((*it).first > idx)    --it;
            acc = (*it).second; //+= idx2len[*it];
            gap_len = (*it).second; //idx2len[*it];
            if (*it != *(anchors.begin()))
                gap_len -= (*(std::prev(it, 1))).second; //idx2len[*(std::prev(it, 1))];
            if (idx >= (*it).first && idx < (*it).first + gap_len)
                return gap::GAP;
        }
        return value_type((*sequence)[idx - acc]);
        */
    }

private:

    constexpr size_type get_gap_length(set_iterator it) const noexcept
    {
        if (it == anchors.begin()) return (*it).second;  //idx2len[*it];
        return (*it).second - (*std::prev(it)).second;  // idx2len[*it] - idx2len[*std::prev(it)];
    }

    // reverse update: increase all anchor gaps AFTER position pos by size, i.e. start position AND size
    void rupdate(size_type const pos, size_type const size)
    {
        if (LOG_LEVEL_AS)
        {
            std::cout << "initial anchor list: \n";
            for (auto anchor : anchors)
                std::cout << "(" << anchor.first << ", " << anchor.second << "), ";
        // post: update succeeding gaps  by shifting position key right, start from right to left to avoid collisions
        // note: >= not supported for std::_Rb_tree_const_iterator
            std::cout << "rupdate(pos=" << pos << ", size=" << pos << ")" << std::endl;
        }
        size_type new_key, new_val;
        //if (LOG_LEVEL_AS) std::cout << "auto it = std::prev(anchors.end()) = " << (*std::prev(anchors.end(), 1)) << std::endl;
        for (auto it = std::prev(anchors.end(), 1); (*it).first > pos;)
        {
            //std::cout << "iteration " << i << std::endl;
            new_key = (*it).first + size;
            new_val = (*it).second + size;  //idx2len[*it] + size;
            if (LOG_LEVEL_AS) std::cout << "reverse update with *it = " << (*it).first << std::endl;
            anchors.insert(gap_t{new_key, new_val});
            //nodeHandler.key() = new_key;
            //idx2len.insert(std::move(nodeHandler));
            //idx2len[new_key] = new_val;
            //anchors.insert(new_key);
            anchors.erase(*it--);
        }
    }

    // forward update: decrease all anchor gaps after position pos by size
    void update(size_type const pos, size_type const size)
    {
        if (LOG_LEVEL_AS) std::cout << "DEBUG: update with pos = " << pos << " and size_del = " << size << std::endl;
        //assert(pos >= size);
        // post: update succeeding gaps  by shifting position key right, start from right to left to avoid collisions
        auto it = anchors.lower_bound(gap_t{pos + size + 1, _});

        while (it != anchors.end())
        {
            if (LOG_LEVEL_AS) {
                std::cout << "it at end: " << (it == anchors.end()) << std::endl;
                std::cout << "\nupdate anchor gap at " << (*it).first << std::endl;
            }
            gap_t gap{(*it).first - size, (*it).second - size};
            anchors.insert(gap);
            if (LOG_LEVEL_AS) std::cout << "erase from anchors: " << (*it).first << std::endl;

            set_iterator it_next = std::next(it);
            anchors.erase(it);
            it = it_next;
            // next line necessary?
            if (it_next == anchors.end()) break;

        }
        if (LOG_LEVEL_AS)
        {
            std::cout << "new anchors: ";
            for (auto a : anchors) std::cout << "(" << a.first << ", " << a.second << "), ";
            std::cout << std::endl;
        }
    }

    // pointer to ungapped, underlying sequence
    inner_type const * sequence{};
    // store virtual gap positions together with the number of gaps until the given position
    std::set<gap_t, gap_compare<gap_t>> anchors{};
    // store accumulated gap lengths corresponding to virtual gap positions
    //std::unordered_map<size_type, size_type> idx2len{};
};

} // namespace seqan3

#endif
