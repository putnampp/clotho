//   Copyright 2015 Patrick Putnam
//
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License.
#ifndef BLOCK_MASKS_HPP__
#define BLOCK_MASKS_HPP__

#include <cassert>

namespace clotho {
namespace utility {

template < class Block >
class block_masks {
public:
    static const unsigned int bits_per_block = sizeof( Block ) * 8;

    typedef Block block_type;
    typedef Block mask_type;

    static const mask_type (&lookupLowBitMask())[bits_per_block] {
        static bool i = init_low_order(low_order_bit_masks);
        assert( i );
        return low_order_bit_masks;
    }

    static const mask_type (&lookupPositionMask())[bits_per_block] {
        static bool i = init_positions(bit_position_masks);
        assert( i );
        return bit_position_masks;
    }

    static mask_type low_order_mask( unsigned int idx ) {
        assert( idx < bits_per_block );
        return ((idx == bits_per_block - 1) ? (mask_type)-1 : (( (mask_type)1 << (idx + 1) ) - 1));
    }

    static mask_type position_mask( unsigned int idx ) {
        assert( idx < bits_per_block );
        return ((mask_type)1 << idx);
    }

private:

    ///  LOW_BIT_MASK => bits_{[0, n]} = 1
    static mask_type low_order_bit_masks[ bits_per_block ];

    /// OFFSET_BIT_MASK => bit_[n] = 1
    static mask_type bit_position_masks[ bits_per_block ];

    static bool init_low_order( mask_type * masks ) {
        for( unsigned int i = 0; i < bits_per_block; ++i ) {
            (*masks++) = low_order_mask( i );
        }
        return true;
    }

    static bool init_positions( mask_type * masks ) {
        for( unsigned int i = 0; i < bits_per_block; ++i ) {
            (*masks++) = position_mask( i );
        }
        return true;
    }
};

template < class Block >
typename block_masks< Block >::mask_type block_masks< Block >::low_order_bit_masks[ block_masks< Block >::bits_per_block ];

template < class Block >
typename block_masks< Block >::mask_type block_masks< Block >::bit_position_masks[ block_masks< Block >::bits_per_block ];
}   // namespace utility {
}   // namespace clotho {
#endif  // BLOCK_MASKS_HPP__
