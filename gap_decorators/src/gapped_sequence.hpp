// wrapper for std::vector<gapped<alphabet_type>> to provide unified interface for simpler benchmark generation

template <typename inner_type>  //=gapped<dna_type>
struct gapped_sequence
{
    using alphabet_type = typename ranges::v3::value_type_t<inner_type>;
    using value_type = inner_type;
    using reference = value_type;
    using const_reference = const reference;
    using iterator = std::vector<inner_type>::iterator;
    using const_iterator = iterator;
    using difference_type = typename ranges::v3::difference_type_t<inner_type>;
    using size_type = typename ranges::v3::size_type_t<inner_type>;

    bool insert_gap(size_type const pos, size_type const size=1)
    {
        vec.insert(vec.begin() + pos, size, inner_type{gap::GAP});
        return true;
    }
    
    bool resize(size_type new_size)
    {
        vec.resize(new_size);
        return true;
    }
    
private:
    std::vector<inner_type> vec;
    
};
