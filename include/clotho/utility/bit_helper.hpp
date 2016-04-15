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
};

template < >
struct bit_helper< 2 > {
    static const unsigned int BITS_PER_BLOCK = 16;
};

template < >
struct bit_helper< 4 > {
    static const unsigned int BITS_PER_BLOCK = 32;
};

template < >
struct bit_helper< 8 > {
    static const unsigned int BITS_PER_BLOCK = 64;
};

template < class T >
struct BitHelper : public bit_helper< sizeof(T) > {
    typedef T block_type;

    static block_type bit_offset( unsigned int idx ) {
        return (block_type) bit_masks[ idx ];
    }

    static block_type low_bit_mask( unsigned int idx ) {
        return (block_type) low_bit_masks[ idx ];
    }
};

}   // namespace utility
}   // namespace clotho

#endif  // BIT_HELPERS_HPP_
