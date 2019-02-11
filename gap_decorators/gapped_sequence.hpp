// wrapper for std::vector<gapped<alphabet_type>> to provide unified interface for simpler benchmark generation
#include <vector>

#include <range/v3/all.hpp>

#include <seqan3/alphabet/gap/gapped.hpp>

namespace seqan3
{

template <typename inner_type> // inner_type requires container concept
struct gapped_sequence
{
    using value_type = typename inner_type::value_type;  //ranges::v3::value_type_t<inner_type>;  // =gapped<dna_type>
    using reference = value_type;
    using const_reference = const reference;
    using iterator = typename inner_type::iterator; //std::vector<value_type>::iterator;
    using const_iterator = iterator;
    using difference_type = typename inner_type::difference_type;  //ranges::v3::difference_type_t<value_type>;
    using size_type = typename inner_type::size_type;  //ranges::v3::size_type_t<value_type>;

    constexpr gapped_sequence()
    {
        _data = inner_type();  //std::vector<value_type>();
    };

    constexpr gapped_sequence(gapped_sequence const &) = default;

    constexpr gapped_sequence & operator=(gapped_sequence const &) = default;

    constexpr gapped_sequence (gapped_sequence && rhs) = default;

    constexpr gapped_sequence & operator=(gapped_sequence && rhs) = default;

    ~gapped_sequence() = default;

    constexpr gapped_sequence(inner_type * sequence)
    {
        for (auto it = sequence->begin(); it < sequence->end(); ++it)
            _data.push_back(*it);
    };

    bool insert_gap(size_type const pos, size_type const size=1)
    {
        _data.insert(_data.begin() + pos, size, value_type{gap::GAP});
        return true;
    }

    bool erase_gap(size_type const pos1, size_type const pos2)
    {
        assert(pos1 < pos2 && pos2 <= this->size());
        // check for consecutive gap
        for (size_type i = pos1; i < pos2; ++i)
            assert(_data.at(i) == gap::GAP);
        _data.erase(_data.begin() + pos1, _data.begin() + pos2);
        return true;
    }

    size_type size()
    {
        return _data.size();
    }

    bool resize(size_type new_size)
    {
        _data.resize(new_size);
        return true;
    }

    constexpr reference operator[](size_type const idx) // const noexcept(noexcept((*host)[pos+n]))
    {
        return _data[idx];
    }

private:
    inner_type _data;

};

}
