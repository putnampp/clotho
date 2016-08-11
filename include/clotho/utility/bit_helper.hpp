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
#ifndef BIT_HELPERS_HPP_
#define BIT_HELPERS_HPP_

#include "clotho/utility/bit_masks.hpp"

namespace clotho {
namespace utility {

template < unsigned char BYTE_COUNT >
struct bit_helper;

template < >
struct bit_helper< 1 > {
    static const unsigned int BITS_PER_BLOCK = 8;
    static const unsigned long long POW2_MOD = LOW_POW2_2;
    static const unsigned char ALL_SET = 0xFF;
    static const unsigned char ALL_UNSET = 0x00;
};

template < >
struct bit_helper< 2 > {
    static const unsigned int BITS_PER_BLOCK = 16;
    static const unsigned long long POW2_MOD = LOW_POW2_3;
    static const unsigned short ALL_SET = 0xFFFF;
    static const unsigned short ALL_UNSET = 0x0000;
};

template < >
struct bit_helper< 4 > {
    static const unsigned int BITS_PER_BLOCK = 32;
    static const unsigned long long POW2_MOD = LOW_POW2_4;
    static const unsigned int ALL_SET = 0xFFFFFFFF;
    static const unsigned int ALL_UNSET = 0x00000000;
};

template < >
struct bit_helper< 8 > {
    static const unsigned int BITS_PER_BLOCK = 64;
    static const unsigned long long POW2_MOD = LOW_POW2_5;
    static const unsigned long long ALL_SET = 0xFFFFFFFFFFFFFFFF;
    static const unsigned long long ALL_UNSET = 0x0000000000000000;
};

template < class T >
struct BitHelper : public bit_helper< sizeof(T) > {
    typedef T block_type;
    typedef bit_helper< sizeof(T) > base_type;
    static const unsigned int   SSE_ALIGNED_SIZE = 128;
    static const unsigned int   BLOCKS_PER_SSE = (SSE_ALIGNED_SIZE / base_type::BITS_PER_BLOCK);

#define MOD_POW2( x ) x & base_type::POW2_MOD

    static block_type bit_offset( size_t idx ) {
        return (block_type) bit_masks[ MOD_POW2( idx ) ];
    }

    static block_type low_bit_mask( size_t idx ) {
        return (block_type) low_bit_masks[ MOD_POW2( idx )  ];
    }

    static size_t   padded_block_count( size_t s ) {
        return (BLOCKS_PER_SSE * ((s / SSE_ALIGNED_SIZE) + ( (s % SSE_ALIGNED_SIZE) ? 1 : 0)));
    }

#undef MOD_POW2
};

}   // namespace utility
}   // namespace clotho

#endif  // BIT_HELPERS_HPP_
