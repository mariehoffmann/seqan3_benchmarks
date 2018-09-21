/*
    Uses an unordered set for the virtual gap positions to retrieve log (k) lookup time.
    Drawback: since the virtual gap positions correspond to the accumulative sum of gap lengths, on update or deletion, the
    worst-case runtime is O(k) for updating subsequent gap positions.
    lower_bound probably does not take advantage of ordering when inserting tuple/pair, because order criterion not obvious.
 */

//#pragma once

#ifndef MYHEADER_H

#define MYHEADER_H

#include <iostream>
#include <set>
#include <stdlib.h>
#include <unordered_map>

#include <range/v3/all.hpp>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>

#define LOG_LEVEL_AS 0

namespace seqan3 {

template <typename inner_type>
struct anchor_set
{

public:

    using alphabet_type = typename ranges::v3::value_type_t<inner_type>;
    using value_type = gapped<alphabet_type>;
    using reference = value_type;
    using size_type = typename ranges::v3::size_type_t<inner_type>;

    constexpr anchor_set()
    {
        data = std::shared_ptr<data_t>(new data_t{});
    };

    constexpr anchor_set(anchor_set const &) = default;

    constexpr anchor_set & operator=(anchor_set const &) = default;

    constexpr anchor_set (anchor_set && rhs) = default;

    constexpr anchor_set & operator=(anchor_set && rhs) = default;

    ~anchor_set() = default;

    constexpr anchor_set(inner_type * sequence): data{new data_t{sequence}} {};


    size_type size() const noexcept
    {

        if (data->anchor_idcs.size())
            return data->idx2len[*(data->anchor_idcs.rbegin())] + data->sequence->size();
        return data->sequence->size();
    }

    // for benchmark only
    bool resize(size_type new_size)
    {
        if (LOG_LEVEL_AS) std::cout << "start resizing ... \n";
        assert(new_size <= this->size());
        //assert(data->sequence->size() > 0);
        for (auto pos = this->size() - 1; pos >= new_size; --pos)
        {
            if (LOG_LEVEL_AS) std::cout << "\tcurrent pos = " << pos << std::endl;
            if ((value_type)(*this)[pos] == gap::GAP)
            {
                if (LOG_LEVEL_AS) std::cout << "\tcase: del gap";
                erase_gap(pos);
            }
            else
            {
                if (LOG_LEVEL_AS) std::cout << "\tcase: del letter";
                data->sequence->resize(data->sequence->size() - 1);
            }
        }
        if (LOG_LEVEL_AS) std::cout << "... final size = " << this->size() << std::endl;
        return true;
    }

    bool insert_gap(size_type const pos, size_type const size=1)
    {
        if (LOG_LEVEL_AS) std::cout << "called insert with pos = " << pos << ", size = " << size << std::endl;
        typename std::set<size_type>::iterator it;
        // case 1: extend previous/surrounding gap
        if ((pos < this->size()) && (((value_type)(*this)[pos] == gap::GAP) or (pos > 0 && (value_type)(*this)[pos-1] == gap::GAP)))
        {
            if (LOG_LEVEL_AS) std::cout << "case: gap extension\n" << std::endl;
            it = data->anchor_idcs.lower_bound(pos);
            if (it == data->anchor_idcs.end() || *it > pos)
            {
                it = std::prev(it);
                if (LOG_LEVEL_AS) std::cout << "case: decrease iterator\n";
            }
            data->idx2len[*it] += size;
            if (LOG_LEVEL_AS)
                std::cout << "updated existing anchor pos = " << *it << " with new acc gap length = " << data->idx2len[*it] << std::endl;
            // merge with successor
            typename std::set<size_type>::iterator it_next = it;
            ++it_next;

            if (*it < *std::prev(data->anchor_idcs.end()) && (*it_next) <= (*it) + size - 1)
            {
                if (LOG_LEVEL_AS) std::cout << "case: merge also with successor (" << *(it_next) << ", " << data->idx2len[*(it_next)] << ") ... \n";
                // extend gap for *it, delete *(it+1)
                data->idx2len[*it] += data->idx2len[*(it_next)];
                data->anchor_idcs.erase(*(it_next));
                data->idx2len.erase(*(it_next));

            }

        }
        // case 2: new anchor gap
        else
        {
            if (LOG_LEVEL_AS) std::cout << "case: new anchor gap: " << std::endl;
            data->anchor_idcs.insert(pos);
            data->idx2len[pos] = size;
            if (LOG_LEVEL_AS) std::cout << "anchor_idx = " << pos << ", len = " << data->idx2len[pos] << std::endl;
            // pre: pos not in anchor set, what's the next lower index?
            it = data->anchor_idcs.find(pos); // return value in case of no lower elem?

            if (LOG_LEVEL_AS)
            {
                for (auto it2 = data->anchor_idcs.begin(); it2 != data->anchor_idcs.end(); ++it2)
                    std::cout << "elem in anchor set: " << (*it2) << std::endl;
            }
            if (LOG_LEVEL_AS) std::cout << "lower elem = " << *it << std::endl;
            // add accumulated gaps from preceeding gap
            if (it != data->anchor_idcs.begin()){
                if (LOG_LEVEL_AS)
                    std::cout << "add acc from previous: " << data->idx2len[*prev(it)] << std::endl;
                data->idx2len[pos] += data->idx2len[*--it];
                if (LOG_LEVEL_AS) std::cout << "accumlated gaps for new gaps: " << data->idx2len[pos] << std::endl;
            }
        }
        // post-processing: reverse update of succeeding gaps
        if (LOG_LEVEL_AS) std::cout << "call rupdate with pos = " << pos << " and size = " << size << std::endl;
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
        /*if ((value_type)(*this)[pos1] != gap::GAP || pos1 >= pos2){
            std::cout << "error no gap at this pos\n";
            std::exit(-1);
            return false;
        }*/
        typename std::set<size_type>::iterator it = data->anchor_idcs.lower_bound(pos1);
        size_type gap_len = get_gap_length(it);

        if (it == data->anchor_idcs.end() || (*it) > pos1) it = std::prev(it);
        if (LOG_LEVEL_AS) std::cout << "lower_bound is: " << *it << std::endl;
        // case 1: complete gap is deleted
        if (LOG_LEVEL_AS) std::cout << "gap_len = " << data->idx2len[pos1] << " and pos2-pos1 = " << pos2-pos1 << std::endl;
        if (((*it) == pos1) && (gap_len == pos2-pos1))
        {
            if (LOG_LEVEL_AS) std::cout << "complete gap del\n";
            data->idx2len.erase(pos1);
            data->anchor_idcs.erase(pos1);
        }
        // case 2: gap to be deleted in tail or larger than 1 (equiv. to shift tail left, i.e. pos remains unchanged)
        else
        {
            if (LOG_LEVEL_AS) std::cout << "gap reduction = " << pos2-pos1 << std::endl;
            data->idx2len[*it] -= pos2-pos1;
            if (LOG_LEVEL_AS) std::cout << "new gap len is " << data->idx2len[*it] << std::endl;
        }

        // post-processing: forward update of succeeding gaps
        if (LOG_LEVEL_AS) std::cout << "call update ...\n";
        update(pos1, pos2-pos1);
        if (LOG_LEVEL_AS) std::cout << "... update done\n";
        return true;
    }

    void set_underlying_sequence(inner_type * sequence) const
    {
        data->sequence = sequence;
    }

    constexpr reference operator[](size_type const idx) const
    {
        assert(idx < size());
        // case 1: no gaps
        if (!data->anchor_idcs.size()) return value_type((*data->sequence)[idx]);
        // case 2: gaps
        typename std::set<size_type>::iterator it = data->anchor_idcs.lower_bound(idx);
        if (data->anchor_idcs.size() && it == data->anchor_idcs.end())  it = std::prev(it);

        size_type acc = 0, gap_len = 0;
        if ((*it) <= idx || (it != data->anchor_idcs.begin() && *(std::prev(it)) <= idx))
        {
            if ((*it) > idx)    --it;
            acc += data->idx2len[*it];
            gap_len = data->idx2len[*it];
            if (*it != *(data->anchor_idcs.begin()))
                gap_len -= data->idx2len[*(std::prev(it, 1))];
            if (idx >= (*it) && idx < (*it) + gap_len)
                return gap::GAP;
        }
        return value_type((*data->sequence)[idx - acc]);
    }

private:

    constexpr size_type get_gap_length(typename std::set<size_type>::iterator it) const noexcept
    {
        if (it == data->anchor_idcs.begin()) return data->idx2len[*it];
        return data->idx2len[*it] - data->idx2len[*std::prev(it)];
    }

    // reverse update: increase all anchor gaps AFTER position pos by size, i.e. start position AND size
    void rupdate(size_type const pos, size_type const size)
    {
        // post: update succeeding gaps  by shifting position key right, start from right to left to avoid collisions
        // note: >= not supported for std::_Rb_tree_const_iterator
        if (LOG_LEVEL_AS) std::cout << "rupdate(pos=" << pos << ", size=" << pos << ")" << std::endl;
        size_type new_key, new_val;
        if (LOG_LEVEL_AS) std::cout << "auto it = std::prev(data->anchor_idcs.end()) = " << (*std::prev(data->anchor_idcs.end(), 1)) << std::endl;
        for (auto it = std::prev(data->anchor_idcs.end(), 1); *it > pos;)
        {
            //std::cout << "iteration " << i << std::endl;
            new_key = (*it) + size;
            new_val = data->idx2len[*it] + size;
            if (LOG_LEVEL_AS) std::cout << "reverse update with *it = " << (*it) << std::endl;
            auto nodeHandler = data->idx2len.extract(*it);

            nodeHandler.key() = new_key;
            if (LOG_LEVEL_AS) {
                std::cout << "\tnew key for anchor " << (*it) << " -> " << nodeHandler.key() << std::endl;
                std::cout << "\tnew val for anchor " << (*it) << " -> " << new_val << std::endl;
            }
            data->idx2len.insert(std::move(nodeHandler));
            data->idx2len[new_key] = new_val;


            if (LOG_LEVEL_AS) std::cout << "\tnew value for anchor " << new_key << ": " << data->idx2len[new_key] + size << " (old key = " << data->idx2len[new_key] << ")" << std::endl;

            //assert(data->idx2len.end() != data->idx2len.find(*it));

            //std::cout << "updated key to accum gap size resolution\n";
            // // todo: best would be a guaranteed in place update, since keys are monotonously increased, no balancing needed
            data->anchor_idcs.insert(new_key);
            data->anchor_idcs.erase(*it--);
            //std::cout << "erased old value, it points to " << (*it) << std::endl;
        }

        if (LOG_LEVEL_AS)
        {
            std::cout << "\tupdated anchors: ";
            for (auto anchor : data->anchor_idcs) std::cout << anchor << ": " << data->idx2len[anchor] << "\t";
            std::cout << std::endl;
            std::cout << "total as size : " << this->size() << std::endl;
        }
    }

    // forward update: decrease all anchor gaps after position pos by size
    void update(size_type const pos, size_type const size)
    {
        assert(pos >= size);
        // post: update succeeding gaps  by shifting position key right, start from right to left to avoid collisions
        auto it = data->anchor_idcs.lower_bound(pos + size + 1);
        if (LOG_LEVEL_AS)
        {
            std::cout << *it << std::endl;
            std::cout << (it == data->anchor_idcs.end()) << std::endl;
            std::cout << (it != data->anchor_idcs.end()) << std::endl;
        }

        if (LOG_LEVEL_AS) {std::cout << "anchors: "; for (auto a : data->anchor_idcs) std::cout << a << ", ";}

        while (it != data->anchor_idcs.end()) //for (int i = 0; it != data->anchor_idcs.end(), i < 5; ++i)
        {
            if (LOG_LEVEL_AS) {
                std::cout << "it at end: " << (it == data->anchor_idcs.end()) << std::endl;
                std::cout << "\nupdate anchor gap at " << *it << std::endl;
            }
            auto nodeHandler = data->idx2len.extract(*it);
            nodeHandler.key() = (*it) - size;

            if (LOG_LEVEL_AS) std::cout << "new hdlr key = " << (*it) - size << std::endl;
            data->idx2len.insert(std::move(nodeHandler));
            data->idx2len[(*it) - size] -= size;
            // // todo: best would be a guaranteed in place update, since keys are monotonously increased, no balancing needed
            if (LOG_LEVEL_AS) std::cout << "insert into anchors: " << (*it) - size << std::endl;
            data->anchor_idcs.insert((*it) - size);
            if (LOG_LEVEL_AS) std::cout << "erase from anchors: " << (*it) << std::endl;

            typename std::set<size_type>::iterator it_next = std::next(it);
            data->anchor_idcs.erase(*it);
            it = it_next;
            if (it_next == data->anchor_idcs.end()) break;

        }
        if (LOG_LEVEL_AS)
        {
            std::cout << "new anchors: ";
            for (auto a : data->anchor_idcs) std::cout << "(" << a << ", " << data->idx2len[a] << "), ";
            std::cout << std::endl;
        }
    }


    struct data_t
    {
        // pointer to ungapped, underlying sequence
        inner_type * sequence{};
        // store virtual gap positions
        std::set<size_type> anchor_idcs{};
        // store accumulated gap lengths corresponding to virtual gap positions
        std::unordered_map<size_type, size_type> idx2len{};

    };
    std::shared_ptr<data_t> data;
};

} // namespace seqan3

#endif
