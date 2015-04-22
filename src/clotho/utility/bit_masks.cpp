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
#include "clotho/utility/bit_masks.hpp"

namespace clotho {
namespace utility {

unsigned int scan_set_bits( unsigned int b, unsigned int * indices ) {
    unsigned int count = 0;
    while( b ) {
        unsigned int tmp = LEAST_SIG_BIT( b );
        indices[ count++ ] = DEBRUIJNBIT_HASH_LOOKUP( tmp );
        b ^= tmp;
    }
    return count;
}

unsigned int scan_set_bits( unsigned long b, unsigned int * indices ) {
    unsigned int count = 0;
    unsigned int lo = (unsigned int) b;
    while( lo ) {
        unsigned int tmp = LEAST_SIG_BIT( lo );
        indices[ count++ ] = DEBRUIJNBIT_HASH_LOOKUP( tmp );
        lo ^= tmp;
    }

    lo = (unsigned int)( b >> 32 );
    while( lo ) {
        unsigned int tmp = LEAST_SIG_BIT( lo );
        indices[ count++ ] = DEBRUIJNBIT_HASH_LOOKUP_HI( tmp );
        lo ^= tmp;
    }
    return count;
}

unsigned int hash_set_bits( unsigned int b, unsigned int * indices ) {
    unsigned int count = 0;
    while( b ) {
        unsigned int tmp = LEAST_SIG_BIT( b );
        indices[ count++ ] = DEBRUIJNBIT_HASH( tmp );
        b ^= tmp;
    }
    return count;
}

unsigned int hash_set_bits( unsigned long b, unsigned int * indices ) {
    unsigned int count = 0;
    unsigned int lo = (unsigned int) b;
    while( lo ) {
        unsigned int tmp = LEAST_SIG_BIT( lo );
        indices[ count++ ] = DEBRUIJNBIT_HASH( tmp );
        lo ^= tmp;
    }

    lo = (unsigned int)( b >> 32 );
    while( lo ) {
        unsigned int tmp = LEAST_SIG_BIT( lo );
        indices[ count++ ] = DEBRUIJNBIT_HASH( tmp ) + 32;
        lo ^= tmp;
    }
    return count;
}

}   // namespace utility {
}   // namespace clotho {
