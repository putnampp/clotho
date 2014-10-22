#ifndef BLOCK_MAP_HPP_
#define BLOCK_MAP_HPP_

#include "clotho/powerset/normalized_key.hpp"

#include "clotho/powerset/next_available_order_tag.hpp"

namespace clotho {
namespace powersets {

template < class Element, class Block, class Tag = next_available_order_tag >
struct block_map {
    typedef Element element_type;
    typedef Block   size_type;

    static const unsigned int       bits_per_block = sizeof( size_type ) * 8;

    size_type operator()( const element_type & elem ) {
        return -1;
    }
};

template < class Element, class Block >
struct block_map< Element, Block, normalized_key< Element > > {
    typedef Element element_type;
    typedef Block   size_type;

    typedef normalized_key< Element > keyer;

    static const unsigned int       bits_per_block = sizeof( size_type ) * 8;
    constexpr static const double   width_per_bin = ((keyer::range) / (double) bits_per_block);

    inline size_type operator()( const element_type & elem ) {
        static keyer k;
        return k(elem) * bits_per_block;
    }
};

}   // namespace powersets
}   // namespace clotho

#endif  // BLOCK_MAP_HPP_
