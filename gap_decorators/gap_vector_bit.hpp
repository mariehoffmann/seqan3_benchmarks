#pragma once

#include <algorithm>
#include <initializer_list>
#include <type_traits>
#include <vector>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
//#include <seqan3/range/container/concept.hpp>
#include <seqan3/range/all.hpp>
#include <seqan3/range/detail/random_access_iterator.hpp>

#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>

#define LOG_LEVEL_GV 0

namespace seqan3 {

    template <typename inner_type>
    //!\cond
    requires alphabet_concept<ranges::v3::value_type_t<inner_type>>
    // see doc, in namespace std and renamed: && random_access_range_concept<inner_type> && sized_range_concept<inner_type>
    //!\endcond
    struct aligned_sequence_adaptor_constant_access
    {

    private:
        //!\privatesection
        using aligned_sequence_t    = aligned_sequence_adaptor_constant_access;
        //using alphabet_t = typename ranges::v3::value_type_t<inner_type>::alphabet_type;
        //!\brief Type of the bit-vector.
        using bit_vector_t          = sdsl::bit_vector<>;
        //!\brief Type of the rank support data structure.
        using rank_1_support_t        = sdsl::rank_support_v5<1>;     // bit_vector_t::rank_1_type;
        //!\brief Type of the 0 select support data structure.
        using select_0_support_t      = sdsl::select_support_mcl<0>;  //bit_vector_t::select_0_type;
        //!\brief Type of the 1 select support data structure.
        // space .2*n, access O(1), alternative: select_support_scan space 64 bits, access O(n)
        using select_1_support_t      = sdsl::select_support_mcl<1>;  //bit_vector_t::select_0_type;
        // TODO: adjust builder function

    public:
        //!\publicsection
        /*!\name Member types
         * \{
         */
        //!\brief Value type of container elements.
        //!\hideinitializer
        using alphabet_type = typename ranges::v3::value_type_t<inner_type>;
        using value_type = gapped<alphabet_type>;

        //!\brief Use reference type defined by container.
        //!\hideinitializer
        using reference = value_type;

        //!\brief Use const reference type provided by container.
        //!\hideinitializer
        using const_reference = const reference;
        // decltype(std::as_const(data_values)

        //!\brief Use random access iterator on container as iterator type.
        //!\hideinitializer
        using iterator = detail::random_access_iterator<aligned_sequence_t>;

        //!\brief Use const random access iterator on container as const iterator type.
        //!\hideinitializer
        using const_iterator = iterator;

        //!\brief Type for distances between iterators is taken from alphabet container.
        //!\hideinitializer
        using difference_type = typename ranges::v3::difference_type_t<inner_type>;

        //!\brief Use alphabet container's size_type as a position.
        //!\hideinitializer
        using size_type = typename ranges::v3::size_type_t<inner_type>;
        //!\}

        /* rule of six */
        /*!\name Constructors, destructor and assignment
         * \{
         */
        // \brief Default constructor.
        constexpr aligned_sequence_adaptor_constant_access()
        {
            data = std::shared_ptr<data_t>(new data_t{});
        };

        //!\brief Default copy constructor.
        constexpr aligned_sequence_adaptor_constant_access(aligned_sequence_adaptor_constant_access const &) = default;

        //!\brief Default copy construction via assignment.
        constexpr aligned_sequence_adaptor_constant_access & operator=(aligned_sequence_adaptor_constant_access const &) = default;

        //!\brief Move constructor.
        constexpr aligned_sequence_adaptor_constant_access (aligned_sequence_adaptor_constant_access && rhs) = default;

        //!\brief Move assignment.
        constexpr aligned_sequence_adaptor_constant_access & operator=(aligned_sequence_adaptor_constant_access && rhs) = default;

        //!\brief Use default deconstructor.
        ~aligned_sequence_adaptor_constant_access() = default;
        //!\}

        //!\brief
        /*!\name Constructors of sequence concept.
         * \{
         */
        //!\brief Construct by single value repeated 'size' times
        //shared_ptr<data_t>
        constexpr aligned_sequence_adaptor_constant_access(inner_type * sequence): data{new data_t{sequence}}
        {
            data->gap_vector = bit_vector_t(data->sequence->size(), 0);
        };
        //!\}

        /*!\name Iterators
         * \{
         */
        //!\brief Return iterator pointing to first element of underlying sequence.
        auto begin() noexcept
        {
            return iterator{*this, 0};
        }

        //!\brief Return iterator pointing to past-the-end element of gapped sequence.
        auto end() noexcept
        {
            return iterator{*this, size()};
        }
        //!\}

        /*!\name Boolean operators
         * \{
         */
        /*!\brief Equality operator for aligned sequences.
         *
         * Two aligned sequences are the same if their literal sequences and gap
         * positions are the same. Note there is no operator== in sdsl-lite and the one
         * in the sdsl master branch compares only addresses, but we are interested in
         * the content. This implementation has best case runtime O(1) and worst-case
         * O(r+s+m) where r,s are the initialization costs of sdsl::rank_1_support_sd,
         * and sdsl::select_support_sd, and m the number of 1 bits.
         */
        bool operator==(aligned_sequence_t & rhs)
        {
            if (data->sequence != rhs.data->sequence || this->size() != rhs.size())
                return false;
            if (data->dirty)
                update_support_structures();
            if (rhs.data->dirty)
                rhs.update_support_structures();
            if (data->rank_1_support.rank(this->size()) != rhs.data->rank_1_support.rank(rhs.size()))
                return false;
            size_type m = data->rank_1_support.rank(this->size());
            for (size_type i = 1; i < m; ++i)
                if (data->select_1_support.select(i) != rhs.data->select_1_support.select(i))
                    return false;
            return true;
        }

        //!\brief Unequality operator for aligned sequences.
        bool operator!=(aligned_sequence_t & rhs)
        {
            return !(*this == rhs);
        }

        //!\brief Swap two aligned sequences and their support structures.
        void swap(aligned_sequence_t & rhs)
        {
            data.swap(rhs.data);
        }

        //!\brief Return gapped sequence length.
        size_type size() noexcept
        {
            // either keep it uptodate, also when resizing underlying sequence (del 0s!) or use rank
            // return data->gap_vector.size();
            if (data->dirty)
                update_support_structures();

            if (LOG_LEVEL_GV)
            {
                std::cout << "enter size ...\n";
            std::cout << "is dirty = " << data->dirty << std::endl;
            std::cout << "data->gap_vector.size() = " << data->gap_vector.size() << std::endl;
            std::cout << "data->rank_1_support.rank(1) = " <<  data->rank_1_support.rank(0) << std::endl;
            }
            return data->rank_1_support.rank(data->gap_vector.size()-1) + data->sequence->size();
        }

        /*!\brief Return the maximal aligned sequence length.
         *
         * The maximal sequence length is limited by either the maximal size of the
         * compressed sequence or the gap vector. Since the sdsl::sd_vector has no
         * max_size() member function, but can be constructed by a bit_vector, we
         * assume the maximal size is restricted by the one of sdsl::bit_vector.
         */
        size_type max_size() const
        {
            return std::min<size_type>(data->sequence->max_size(), sdsl::bit_vector{}.max_size());
        }

        //!\brief An aligned sequence is empty if it contains no alphabet letters or gaps.
        bool empty() const
        {
            return data->sequence->empty() && data->gap_vector.size() == 0;
        }
        //!\}

        /*!\name Sequence concept support.
         * \{
         */
        /*!\brief Insert a single value at a given iterator position.
         *
         * Elements right of the newly inserted elemented are shifted right. The
         * returned iterator points to the position of the insertion.
         */

        /*!\brief Insert a value multiple times at an iterator position.
         *
         * The returned iterator points to the position of the left-most inserted
         * element.
         */
        iterator insert_gap(iterator it, size_type size=1)
        {
            size_type pos = static_cast<size_type>(it - detail::random_access_iterator<aligned_sequence_t>());
            assert(insert_gap(pos, size));
            return it;
        }

        /*!\brief Insert value at a position 'size' times.
         *
         * Return false if given position exceeds the current size by one.
         */
        bool insert_gap(size_type const pos, size_type const size=1)
        {
            if (LOG_LEVEL_GV) std::cout << "enter insert_gap with pos = " << pos << ", and size = " << size << std::endl;
            if (pos > this->size())
                return false;
            if (data->dirty)
                update_support_structures();
            // rank queries on empty sd_vector throws assertion
            size_type m = (!this->size()) ? size : data->rank_1_support.rank(this->size()) + size;
            if (LOG_LEVEL_GV) std::cout << "init builder with new amount of set bits: " << m << std::endl;
            data->bit_vector.resize(this->size() + size);

            // shift suffix, note that we need i to be a signed integer for the case
            // that the aligned sequence was empty
            if (LOG_LEVEL_GV) std::cout << "iterate over suffix from i = " << this->size() - 1 << " to " << static_cast<signed int>(pos) << std::endl;
            for (size_t i = this->size()-1; i >= pos; --i)
                data->gap_vector.set_int(i+size, data->gap_vector[i]);

            // insert gap
            for (size_type i = pos; i < pos+size; ++i)
                data->gap_vector.set_int(i, 1);

            if (LOG_LEVEL_GV)
            {
                std::cout << "gap vector after insertion: ";
                for (auto bit : data->gap_vector) std::cout << bit;
                std::cout << std::endl;
            }
            return true;
        }

        /*\brief Erase element at the iterator's position.
         *
         * Return iterator to past-the-erased element.
         */
        iterator erase_gap(iterator const it)
        {
            assert(erase_gap(static_cast<size_type>(it - iterator{})));
            return it;
        }

        /*!\brief Erase element at a given index.
         *
         * Return false if index exceeds current size minus one.
         * Worst-case runtime is O(r+s+m). r,s re-initialization of rank and select
         * support, m number of set bits.
         */
        bool erase_gap(size_type const pos)
        {
            if (!data->gap_vector[pos])
                return false;
            return erase_gap(pos, pos+1);
        }

        //!\brief Erase all gaps in range pos1 and pos2 (exclusive). Gaps
        // right-hand of 2nd iterator are shifted by the number of gaps deleted.
        bool erase_gap(size_type const pos1, size_type const pos2)
        {
            assert(pos1 <= pos2);
            //std::cout << "enter erase_gap with pos1 = " << pos1 << ", pos2 = " << pos2 << std::endl;
            if (pos1 >= size() || pos2 > size())
                return false;
            if (data->dirty)
                update_support_structures();

            // shift suffix at it2 by the size_type m_del it2-it1
            size_type pos_del = pos2 - pos1;
            for (size_type i = pos1; i < std::min<size_type>(this->size(), pos2); ++i)
                data->gap_vector.set_int(i, data->gap_vector[i - pos_del]);

            data->gap_vector.resize(this->size() - pos_del);
            data->dirty = true;
            return true;
        }

        //!\brief Erase all gaps falling into range given by first and second
        // iterator (exclusive). Gaps right-hand of 2nd iterator are shifted by
        // the number of gaps deleted.
        iterator erase_gap(iterator const it1, iterator const it2)
        {
            size_type pos1 = static_cast<size_type>(it1 - iterator{});
            size_type pos2 = static_cast<size_type>(it2 - iterator{});
            assert(erase_gap(pos1, pos2));
            return it1;
        }

        //!\brief Append gap symbol to the aligned sequence.
        // Note: there is no resize for sd_vector, it has to be re-initialized with
        // sd_vector_builder.
        void push_back()
        {
            assert(max_size() >= size() + 1);
            data->gap_vector.resize(this->size() + 1, 1);
            data->dirty = true;
        }

        //!\brief Delete last gap symbol if set, else return false.
        //
        // Worst-case runtime O(r+s+m), r, s rank and select support initialization,
        // m number of 1 bits.
        // TODO: behaviour - pop last gap or (wherever it occurs) or only query
        // last position and remove if gap?
        bool pop_back()
        {
            assert(this->size() > 0);
            if (!data->gap_vector[size() - 1])
                return false;
            data->gap_vector.resize(this->size() - 1);
            data->dirty = true
            return true;
        }

        //!\brief Clear gaps in bit vector. Alphabet sequence remains unchanged.
        void clear()
        {
            data->gap_vector.assign(data->sequence->size(), 0);
            data->dirty = true;
        }

        //!\brief Return first symbol of aligned sequence.
        reference front()
        {
            assert(size() > 0u);
            return (*this)[0];
        }

        //!\brief Return last symbol of aligned sequence.
        reference back()
        {
            assert(size() > 0u);
            return (*this)[size()-1];
        }
        //!\}

        /*!\name Sequence getter and setter.
         * \{
         */
        //!\brief Return pointer to gap-free sequence.
        inner_type * get_underlying_sequence() const
        {
            return data->sequence;
        }

        //!\brief Set pointer to ungapped sequence.
        void set_underlying_sequence(inner_type * sequence) const
        {
            data->sequence = sequence;
            data->gap_vector.assign(sequence->size(), 0);
            data->dirty = true;
        }
        //!\}

        /*\brief Map a compressed representation index to the aligned sequence index.
         *
         * Note that a i-th 0 corresponds to the i-th position in the compressed sequence.
         * We therefore can directly use the select support structure to map to the
         * aligned sequence space.
         * E.g. '--TA-TA--' with input pos = 2 returns select<0>(2) = 5.
         */
        size_type map_to_aligned_position(size_type const idx)
        {
            if (idx >= size())
                throw std::out_of_range{"Trying to access element behind the last in aligned_sequence."};
            if (data->dirty)
                update_support_structures();
            return data->select_0_support.select(idx + 1);
        }

        /*!\brief Map from gapped sequence position to index of compressed
         * sequence representation.
         *
         * We use the rank support structure indexed for counting 1s and substract
         * the number gaps in [0; position_gap] from the input position.
         * E.g.           aligned sequence  | - A - - T
         *                 position_gapped  | 0 1 2 3 4
         *      map_to_underlying_position  |-1 0 0 0 1
         * Note that a gap is mapped to the same position than the next preceeding
         * non-gap symbol.
         */
        difference_type map_to_underlying_position(size_type const position_gapped)
        {
            if (data->dirty)
                update_support_structures();
            return static_cast<difference_type>(position_gapped) -
            static_cast<difference_type>(data->rank_1_support.rank(std::min<size_type>(position_gapped+1, this->size())));
        }
        //!\}

        /*!\name Random access sequence concept support.
         * \{
         */
        //!\brief Return reference to aligned sequence for given index.
        constexpr reference operator[](size_type const idx) // const noexcept(noexcept((*host)[pos+n]))
        {
            assert(idx < size());
            if (!data->gap_vector[idx]){
                size_type const pos = map_to_underlying_position(idx);
                return value_type((*data->sequence)[pos]);
            }
            return gap::GAP;
        }

        //!\brief Return reference to aligned sequence for given index.
        value_type at(size_type const idx)
        {
            if (idx >= size())
                throw std::out_of_range{"Trying to access element behind the last in aligned_sequence."};
            return (value_type)(*this)[idx];
        }
        //!\}

        bool resize(size_type new_size)
        {
            if (LOG_LEVEL_GV) std::cout << "start resizing ... ";
            assert(new_size <= this->size());
            //assert(data->sequence->size() > 0);
            for (auto pos = this->size() - 1; pos >= new_size; --pos)
            {
                if (LOG_LEVEL_GV) std::cout << "current pos = " << pos << ", gseq[pos] = " << (value_type)(*this)[pos] << std::endl;
                if ((value_type)(*this)[pos] == gap::GAP)
                {
                    erase_gap(pos);
                }
                else
                {
                    data->sequence->resize(data->sequence->size() - 1);
                }
            }
            data->gap_vector.resize(new_size);
            if (LOG_LEVEL_GV) std::cout << "... final size = " << this->size() << std::endl;
            return true;
        }

    private:
        //!\privatesection
        //!\brief Structure for storing a sequence and gap information and helper functions.

        struct data_t
        {
            /*!\brief Pointer to where the ungapped sequence is stored.
             *
             * The ungapped sequence is the original sequence of an ungapped alphabet type.
             * If the alphabet type allows gap symbols, these are treated as normal symbols.
             * Only gaps inserted via this interface are stored in a bit vector.
             * Per default it is a null pointer.
             */
            inner_type * sequence{};

            /*!\brief Where the gapped sequence is stored.
             *
             * Gap presentation of the aligned sequence (1: gap, 0: alphabet letter).
             * The total length is sequence size + number of gaps and therefore corresponds
             * to the aligned sequence size.
             */
            bit_vector_t gap_vector = bit_vector_t();

            /*!\brief Rank support structure for projection into ungapped sequence space.
             *
             * The rank of position i is number 1s in the prefix 0..i-1. This corresponds
             * to the number of gaps in the aligned sequence. Rank queries are answered in
             * constant time.
             */
            //!\hideinitializer
            rank_1_support_t rank_1_support{};
            /*!\brief Select support structure for projection into gap space.
             *
             * Select i returns the gap vector position of the i-th zero in constant time.
             */
            //!\hideinitializer
            select_0_support_t select_0_support{};

            /*!\brief Select support structure for copying set bits in O(m).
             *
             * Select i returns the gap vector position of the i-th 1 in constant time.
             */
            //!\hideinitializer
            select_1_support_t select_1_support{};

            //!\brief Flag to indicate whether sd_vector support structures needs to be
            // updated before executing rank or select queries.
            bool dirty{true};
        };

        std::shared_ptr<data_t> data;

        //!\brief Re-initialize rank and select support structures of bit_vector.
        void update_support_structures()
        {
            data->rank_1_support = sdsl::rank_support_v5<1>(&data->gap_vector);
            data->select_0_support = sdsl::select_support_mcl<0>(&data->gap_vector);
            data->select_1_support = sdsl::select_support_mcl<1>(&data->gap_vector);
            data->dirty = false;
        }
    };

    //!\brief Global swap function.
    template <typename inner_type, char gap_symbol = '_'>
    void swap (aligned_sequence_adaptor_constant_access<inner_type> & lhs, aligned_sequence_adaptor_constant_access<inner_type> & rhs)
    {
        lhs.swap(rhs);
    }

} // namespace seqan3
