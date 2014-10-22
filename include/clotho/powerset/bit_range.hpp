#ifndef CLOTHO_BIT_RANGE_HPP_
#define CLOTHO_BIT_RANGE_HPP_

#include "clotho/powerset/key_range.hpp"

namespace clotho {
namespace powersets {

template < class Block, class Range >
struct bit_range : public key_range < Range > {
    typedef Block block_type;

    static const unsigned int       bits_per_block = sizeof( block_type ) * 8;
    constexpr static const double   width_per_bin = ((key_range< Range >::range) / (double) bits_per_block);
};

}   // namespace powersets {
}   // namespace clotho {

#endif  // CLOTHO_BIT_RANGE_HPP_
