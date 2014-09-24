#ifndef BLOCK_MASKS_HPP_
#define BLOCK_MASKS_HPP_

namespace clotho {
namespace utility {

/**
 * Block masks
 *
 * Object with lookup tables for common masks used in bit-based algorithms.
 *
 */
template < class Block, unsigned int N = (sizeof( Block ) * 8) >
class block_masks  {
public:
    static const unsigned int bits_per_block = N;

    typedef Block block_type;
    typedef Block mask_type;

    ///  LOW_BIT_MASK => bits_{[0, n]} = 1
    static const mask_type low_order_bit_masks[ bits_per_block ];

    /// OFFSET_BIT_MASK => bit_[n] = 1
    static const mask_type bit_position_masks[ bits_per_block ];

    static mask_type low_order_mask( unsigned int idx ) {
        return (( 1 << (idx + 1) ) - 1);
    }

    static mask_type position_mask( unsigned int idx ) {
        return (1 << idx);
    }
};

/// OFFSET_BIT_MASK => bit_[n] = 1
#define OFFSET_BIT_MASK( n ) (1UL << n)

///  LOW_BIT_MASK => bits_{[0, n]} = 1
#define LOW_BIT_MASK( n ) ((1UL << (n + 1) ) - 1)

#define SPEC64 block_masks< Block, 64 >
#define SPEC32 block_masks< Block, 32 >
#define SPEC16 block_masks< Block, 16 >
#define SPEC8 block_masks< Block, 8 >

template < class Block >
const typename SPEC64::mask_type  SPEC64::low_order_bit_masks[ SPEC64::bits_per_block ] = {
     LOW_BIT_MASK( 0 ), LOW_BIT_MASK( 1 ), LOW_BIT_MASK( 2 ), LOW_BIT_MASK( 3 ),
     LOW_BIT_MASK( 4 ), LOW_BIT_MASK( 5 ), LOW_BIT_MASK( 6 ), LOW_BIT_MASK( 7 ),
     LOW_BIT_MASK( 8 ), LOW_BIT_MASK( 9 ), LOW_BIT_MASK( 10 ), LOW_BIT_MASK( 11 ),
     LOW_BIT_MASK( 12 ), LOW_BIT_MASK( 13 ), LOW_BIT_MASK( 14 ), LOW_BIT_MASK( 15 ),
     LOW_BIT_MASK( 16 ), LOW_BIT_MASK( 17 ), LOW_BIT_MASK( 18 ), LOW_BIT_MASK( 19 ),
     LOW_BIT_MASK( 20 ), LOW_BIT_MASK( 21 ), LOW_BIT_MASK( 22 ), LOW_BIT_MASK( 23 ),
     LOW_BIT_MASK( 24 ), LOW_BIT_MASK( 25 ), LOW_BIT_MASK( 26 ), LOW_BIT_MASK( 27 ),
     LOW_BIT_MASK( 28 ), LOW_BIT_MASK( 29 ), LOW_BIT_MASK( 30 ), LOW_BIT_MASK( 31 ),
     LOW_BIT_MASK( 32 ), LOW_BIT_MASK( 33 ), LOW_BIT_MASK( 34 ), LOW_BIT_MASK( 35 ),
     LOW_BIT_MASK( 36 ), LOW_BIT_MASK( 37 ), LOW_BIT_MASK( 38 ), LOW_BIT_MASK( 39 ),
     LOW_BIT_MASK( 40 ), LOW_BIT_MASK( 41 ), LOW_BIT_MASK( 42 ), LOW_BIT_MASK( 43 ),
     LOW_BIT_MASK( 44 ), LOW_BIT_MASK( 45 ), LOW_BIT_MASK( 46 ), LOW_BIT_MASK( 47 ),
     LOW_BIT_MASK( 48 ), LOW_BIT_MASK( 49 ), LOW_BIT_MASK( 50 ), LOW_BIT_MASK( 51 ),
     LOW_BIT_MASK( 52 ), LOW_BIT_MASK( 53 ), LOW_BIT_MASK( 54 ), LOW_BIT_MASK( 55 ),
     LOW_BIT_MASK( 56 ), LOW_BIT_MASK( 57 ), LOW_BIT_MASK( 58 ), LOW_BIT_MASK( 59 ),
     LOW_BIT_MASK( 60 ), LOW_BIT_MASK( 61 ), LOW_BIT_MASK( 62 ), 0xFFFFFFFFFFFFFFFF 
};

template < class Block >
const typename SPEC64::mask_type  SPEC64::bit_position_masks[ SPEC64::bits_per_block ] = {
    OFFSET_BIT_MASK( 0 ), OFFSET_BIT_MASK( 1 ), OFFSET_BIT_MASK( 2 ), OFFSET_BIT_MASK( 3 ),
    OFFSET_BIT_MASK( 4 ), OFFSET_BIT_MASK( 5 ), OFFSET_BIT_MASK( 6 ), OFFSET_BIT_MASK( 7 ),
    OFFSET_BIT_MASK( 8 ), OFFSET_BIT_MASK( 9 ), OFFSET_BIT_MASK( 10 ), OFFSET_BIT_MASK( 11 ),
    OFFSET_BIT_MASK( 12 ), OFFSET_BIT_MASK( 13 ), OFFSET_BIT_MASK( 14 ), OFFSET_BIT_MASK( 15 ),
    OFFSET_BIT_MASK( 16 ), OFFSET_BIT_MASK( 17 ), OFFSET_BIT_MASK( 18 ), OFFSET_BIT_MASK( 19 ),
    OFFSET_BIT_MASK( 20 ), OFFSET_BIT_MASK( 21 ), OFFSET_BIT_MASK( 22 ), OFFSET_BIT_MASK( 23 ),
    OFFSET_BIT_MASK( 24 ), OFFSET_BIT_MASK( 25 ), OFFSET_BIT_MASK( 26 ), OFFSET_BIT_MASK( 27 ),
    OFFSET_BIT_MASK( 28 ), OFFSET_BIT_MASK( 29 ), OFFSET_BIT_MASK( 30 ), OFFSET_BIT_MASK( 31 ),
    OFFSET_BIT_MASK( 32 ), OFFSET_BIT_MASK( 33 ), OFFSET_BIT_MASK( 34 ), OFFSET_BIT_MASK( 35 ),
    OFFSET_BIT_MASK( 36 ), OFFSET_BIT_MASK( 37 ), OFFSET_BIT_MASK( 38 ), OFFSET_BIT_MASK( 39 ),
    OFFSET_BIT_MASK( 40 ), OFFSET_BIT_MASK( 41 ), OFFSET_BIT_MASK( 42 ), OFFSET_BIT_MASK( 43 ),
    OFFSET_BIT_MASK( 44 ), OFFSET_BIT_MASK( 45 ), OFFSET_BIT_MASK( 46 ), OFFSET_BIT_MASK( 47 ),
    OFFSET_BIT_MASK( 48 ), OFFSET_BIT_MASK( 49 ), OFFSET_BIT_MASK( 50 ), OFFSET_BIT_MASK( 51 ),
    OFFSET_BIT_MASK( 52 ), OFFSET_BIT_MASK( 53 ), OFFSET_BIT_MASK( 54 ), OFFSET_BIT_MASK( 55 ),
    OFFSET_BIT_MASK( 56 ), OFFSET_BIT_MASK( 57 ), OFFSET_BIT_MASK( 58 ), OFFSET_BIT_MASK( 59 ),
    OFFSET_BIT_MASK( 60 ), OFFSET_BIT_MASK( 61 ), OFFSET_BIT_MASK( 62 ), OFFSET_BIT_MASK( 63 )
};

template < class Block >
const typename SPEC32::mask_type  SPEC32::low_order_bit_masks[ SPEC32::bits_per_block ] = {
     LOW_BIT_MASK( 0 ), LOW_BIT_MASK( 1 ), LOW_BIT_MASK( 2 ), LOW_BIT_MASK( 3 ),
     LOW_BIT_MASK( 4 ), LOW_BIT_MASK( 5 ), LOW_BIT_MASK( 6 ), LOW_BIT_MASK( 7 ),
     LOW_BIT_MASK( 8 ), LOW_BIT_MASK( 9 ), LOW_BIT_MASK( 10 ), LOW_BIT_MASK( 11 ),
     LOW_BIT_MASK( 12 ), LOW_BIT_MASK( 13 ), LOW_BIT_MASK( 14 ), LOW_BIT_MASK( 15 ),
     LOW_BIT_MASK( 16 ), LOW_BIT_MASK( 17 ), LOW_BIT_MASK( 18 ), LOW_BIT_MASK( 19 ),
     LOW_BIT_MASK( 20 ), LOW_BIT_MASK( 21 ), LOW_BIT_MASK( 22 ), LOW_BIT_MASK( 23 ),
     LOW_BIT_MASK( 24 ), LOW_BIT_MASK( 25 ), LOW_BIT_MASK( 26 ), LOW_BIT_MASK( 27 ),
     LOW_BIT_MASK( 28 ), LOW_BIT_MASK( 29 ), LOW_BIT_MASK( 30 ), 0xFFFFFFFF 
};

template < class Block >
const typename SPEC32::mask_type  SPEC32::bit_position_masks[ SPEC32::bits_per_block ] = {
    OFFSET_BIT_MASK( 0 ), OFFSET_BIT_MASK( 1 ), OFFSET_BIT_MASK( 2 ), OFFSET_BIT_MASK( 3 ),
    OFFSET_BIT_MASK( 4 ), OFFSET_BIT_MASK( 5 ), OFFSET_BIT_MASK( 6 ), OFFSET_BIT_MASK( 7 ),
    OFFSET_BIT_MASK( 8 ), OFFSET_BIT_MASK( 9 ), OFFSET_BIT_MASK( 10 ), OFFSET_BIT_MASK( 11 ),
    OFFSET_BIT_MASK( 12 ), OFFSET_BIT_MASK( 13 ), OFFSET_BIT_MASK( 14 ), OFFSET_BIT_MASK( 15 ),
    OFFSET_BIT_MASK( 16 ), OFFSET_BIT_MASK( 17 ), OFFSET_BIT_MASK( 18 ), OFFSET_BIT_MASK( 19 ),
    OFFSET_BIT_MASK( 20 ), OFFSET_BIT_MASK( 21 ), OFFSET_BIT_MASK( 22 ), OFFSET_BIT_MASK( 23 ),
    OFFSET_BIT_MASK( 24 ), OFFSET_BIT_MASK( 25 ), OFFSET_BIT_MASK( 26 ), OFFSET_BIT_MASK( 27 ),
    OFFSET_BIT_MASK( 28 ), OFFSET_BIT_MASK( 29 ), OFFSET_BIT_MASK( 30 ), OFFSET_BIT_MASK( 31 )
};

template < class Block >
const typename SPEC16::mask_type  SPEC16::low_order_bit_masks[ SPEC16::bits_per_block ] = {
     LOW_BIT_MASK( 0 ), LOW_BIT_MASK( 1 ), LOW_BIT_MASK( 2 ), LOW_BIT_MASK( 3 ),
     LOW_BIT_MASK( 4 ), LOW_BIT_MASK( 5 ), LOW_BIT_MASK( 6 ), LOW_BIT_MASK( 7 ),
     LOW_BIT_MASK( 8 ), LOW_BIT_MASK( 9 ), LOW_BIT_MASK( 10 ), LOW_BIT_MASK( 11 ),
     LOW_BIT_MASK( 12 ), LOW_BIT_MASK( 13 ), LOW_BIT_MASK( 14 ), 0xFFFF 
};

template < class Block >
const typename SPEC16::mask_type  SPEC16::bit_position_masks[ SPEC16::bits_per_block ] = {
    OFFSET_BIT_MASK( 0 ), OFFSET_BIT_MASK( 1 ), OFFSET_BIT_MASK( 2 ), OFFSET_BIT_MASK( 3 ),
    OFFSET_BIT_MASK( 4 ), OFFSET_BIT_MASK( 5 ), OFFSET_BIT_MASK( 6 ), OFFSET_BIT_MASK( 7 ),
    OFFSET_BIT_MASK( 8 ), OFFSET_BIT_MASK( 9 ), OFFSET_BIT_MASK( 10 ), OFFSET_BIT_MASK( 11 ),
    OFFSET_BIT_MASK( 12 ), OFFSET_BIT_MASK( 13 ), OFFSET_BIT_MASK( 14 ), OFFSET_BIT_MASK( 15 )
};

template < class Block >
const typename SPEC8::mask_type  SPEC8::low_order_bit_masks[ SPEC8::bits_per_block ] = {
     LOW_BIT_MASK( 0 ), LOW_BIT_MASK( 1 ), LOW_BIT_MASK( 2 ), LOW_BIT_MASK( 3 ),
     LOW_BIT_MASK( 4 ), LOW_BIT_MASK( 5 ), LOW_BIT_MASK( 6 ), 0xFF
};

template < class Block >
const typename SPEC8::mask_type  SPEC8::bit_position_masks[ SPEC8::bits_per_block ] = {
    OFFSET_BIT_MASK( 0 ), OFFSET_BIT_MASK( 1 ), OFFSET_BIT_MASK( 2 ), OFFSET_BIT_MASK( 3 ),
    OFFSET_BIT_MASK( 4 ), OFFSET_BIT_MASK( 5 ), OFFSET_BIT_MASK( 6 ), OFFSET_BIT_MASK( 7 )
};

#undef SPEC64
#undef SPEC32
#undef SPEC16
#undef SPEC8


}    namespace utility
}    namespace clotho

#endif  // BLOCK_MASKS_HPP_
