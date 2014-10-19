#ifndef BLOCK_MAP_HPP_
#define BLOCK_MAP_HPP_


namespace clotho {
namespace powersets {

template < class Element, class Block >
struct block_map {
    typedef Element element_type;
    typedef Block   size_type;

    static const unsigned int       bits_per_block = sizeof( size_type ) * 8;

    static constexpr double         min_bin_value = 0.0;
    static constexpr double         max_bin_value = 1.0;

    constexpr static const double   width_per_bin = ((max_bin_value - min_bin_value) / (double) bits_per_block);

    size_type operator()( const element_type & elem ) {
        return 0;
    }
};

}   // namespace powersets
}   // namespace clotho

#endif  // BLOCK_MAP_HPP_
