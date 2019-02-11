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

#define LOG_LEVEL_AS2 0
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

    explicit constexpr anchor_set(inner_type const & sequence): sequence{&sequence} {};

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

    // for benchmark only
    /*
    bool resize(size_type new_size)
    {
        if (LOG_LEVEL_AS2)
            std::cout << "enter resize with new_size = " << new_size << std::endl;
        assert(new_size <= this->size());
        //assert(sequence->size() > 0);
        if (LOG_LEVEL_AS2)
        {
            std::cout << "initial aseq: ";
            for (size_type i = 0; i < this->size(); ++i)
                std::cout << (value_type)(*this)[i];
            std::cout << std::endl;
            std::cout << "initial seq len: " << sequence->size() << std::endl;
        }
        for (auto pos = this->size() - 1; pos >= new_size; --pos)
        {
            if (LOG_LEVEL_AS2) std::cout << "query as[" << pos << "]\n";
            if ((value_type)(*this)[pos] == gap::GAP)
            {
                if (LOG_LEVEL_AS2) std::cout << "gap case\n";
                erase_gap(pos);
            }
            else
            {
                if (LOG_LEVEL_AS2) std::cout << "letter case, resize sequence with " << sequence->size()-1 << std::endl;
                sequence->resize(sequence->size() - 1);
            }
        }
        return true;
    }
    */

    bool insert_gap(size_type const pos, size_type const size=1)
    {
        if (LOG_LEVEL_AS2) std::cout << "called insert with pos = " << pos << ", size = " << size << std::endl;
        //typename std::set<gap_t, gap_compare<gap_t>>::iterator it = anchors.begin();
        set_iterator it = anchors.begin();
        //auto it = anchors.begin();
        // case 1: extend previous/surrounding gap, 'or' instead of '||' never threw an error!?
        if ((pos < this->size()) && (((value_type)(*this)[pos] == gap::GAP) || (pos > 0 && (value_type)(*this)[pos-1] == gap::GAP)))
        {
            if (LOG_LEVEL_AS2) std::cout << "case: gap extension\n" << std::endl;
            it = anchors.lower_bound(gap_t{pos, _});
            if (it == anchors.end() || (*it).first > pos)
            {
                it = std::prev(it);
                if (LOG_LEVEL_AS2) std::cout << "case: decrease iterator\n";
            }
            gap_t gap{(*it).first, (*it).second + size};
//            idx2len[*it] += size;
            if (LOG_LEVEL_AS2)
                std::cout << "updated existing anchor pos = " << (*it).first << " with new acc gap length = " << (*it).second << std::endl;
            // merge with successor
            //typename std::set<gap_t, gap_compare>::iterator it_next = it;
            auto it_next = it;
            ++it_next;

            // How does this happen?
            if ((*it) < (*std::prev(anchors.end())) && (*it_next).first <= (*it).first + size - 1)
            {
                if (LOG_LEVEL_AS2) std::cout << "STATUS: Should not go HERE!\n";
                if (LOG_LEVEL_AS2) std::cout << "case: merge also with successor (" << (*it_next).first << ", " << (*it_next).second << ") ... \n";
                // extend gap for *it, delete *(it+1)
                //gap_t gap{(*it).first, (*it).second + (*it_next).second};
                gap.second += (*it_next).second;
                anchors.erase(it_next);
                /*idx2len[*it] += idx2len[*(it_next)];
                anchors.erase(*(it_next));
                idx2len.erase(*(it_next));*/
            }
            anchors.erase(it);
            anchors.insert(gap);

        }
        // case 2: new anchor gap
        else
        {
            if (LOG_LEVEL_AS2) std::cout << "case: new anchor gap: " << std::endl;
            gap_t gap{pos, size};
//            idx2len[pos] = size;
            // pre: pos not in anchor set, what's the next lower index?
            it = anchors.find(gap_t{pos, _}); // return value in case of no lower elem?
            // add accumulated gaps from preceeding gap
            // TODO: same for anchor_set.hpp
            if (it != anchors.begin() && it != anchors.end()){
                gap.second += (*--it).second;
                //idx2len[pos] += idx2len[*--it];
            }
            anchors.insert(gap);
        }
        // post-processing: reverse update of succeeding gaps
        if (LOG_LEVEL_AS2) std::cout << "call rupdate with pos = " << pos << " and size = " << size << std::endl;
        rupdate(pos, size);
        return true;
    }


    bool erase_gap(size_type const pos)
    {
        return erase_gap(pos, pos+1);
    }

    bool erase_gap(size_type const pos1, size_type const pos2)
    {

        if (LOG_LEVEL_AS2) std::cout << "call erase_gap at pos1 = " << pos1 << ", pos2 = " << pos2 << std::endl;
        /*if ((value_type)(*this)[pos1] != gap::GAP || pos1 >= pos2){
            std::cout << "error no gap at this pos\n";
            std::exit(-1);
            return false;
        }*/
        set_iterator it = anchors.lower_bound(gap_t{pos1, _});
        size_type gap_len = get_gap_length(it);

        if (it == anchors.end() || (*it).first > pos1)
            it = std::prev(it);

        // check for contiguous gap
        assert((*it).first <= pos1 && (*it).first + gap_len >= pos2);
        
        // case 1: complete gap is deleted
        if (((*it).first == pos1) && (gap_len == pos2-pos1))
        {
            if (LOG_LEVEL_AS2) std::cout << "complete gap del\n";
            anchors.erase(it);
        }
        // case 2: gap to be deleted in tail or larger than 1 (equiv. to shift tail left, i.e. pos remains unchanged)
        else
        {
            if (LOG_LEVEL_AS2) std::cout << "gap reduction = " << pos2-pos1 << std::endl;
            //idx2len[*it] -= pos2-pos1;
            gap_t gap{(*it).first, (*it).second - pos2 + pos1};
            // TODO: emplace better?
            anchors.erase(it);
            anchors.insert(gap);
        }

        // post-processing: forward update of succeeding gaps
        if (LOG_LEVEL_AS2) std::cout << "call update ...\n";
        update(pos1, pos2-pos1);
        if (LOG_LEVEL_AS2) std::cout << "... update done\n";
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
        if (LOG_LEVEL_AS2)
        {
            std::cout << "initial anchor list: \n";
            for (auto anchor : anchors)
                std::cout << "(" << anchor.first << ", " << anchor.second << "), ";
        // post: update succeeding gaps  by shifting position key right, start from right to left to avoid collisions
        // note: >= not supported for std::_Rb_tree_const_iterator
            std::cout << "rupdate(pos=" << pos << ", size=" << pos << ")" << std::endl;
        }
        size_type new_key, new_val;
        //if (LOG_LEVEL_AS2) std::cout << "auto it = std::prev(anchors.end()) = " << (*std::prev(anchors.end(), 1)) << std::endl;
        for (auto it = std::prev(anchors.end(), 1); (*it).first > pos;)
        {
            //std::cout << "iteration " << i << std::endl;
            new_key = (*it).first + size;
            new_val = (*it).second + size;  //idx2len[*it] + size;
            if (LOG_LEVEL_AS2) std::cout << "reverse update with *it = " << (*it).first << std::endl;
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
        if (LOG_LEVEL_AS2) std::cout << "DEBUG: update with pos = " << pos << " and size_del = " << size << std::endl;
        //assert(pos >= size);
        // post: update succeeding gaps  by shifting position key right, start from right to left to avoid collisions
        auto it = anchors.lower_bound(gap_t{pos + size + 1, _});

        while (it != anchors.end())
        {
            if (LOG_LEVEL_AS2) {
                std::cout << "it at end: " << (it == anchors.end()) << std::endl;
                std::cout << "\nupdate anchor gap at " << (*it).first << std::endl;
            }
            gap_t gap{(*it).first - size, (*it).second - size};
            anchors.insert(gap);
            if (LOG_LEVEL_AS2) std::cout << "erase from anchors: " << (*it).first << std::endl;

            set_iterator it_next = std::next(it);
            anchors.erase(it);
            it = it_next;
            // next line necessary?
            if (it_next == anchors.end()) break;

        }
        if (LOG_LEVEL_AS2)
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
