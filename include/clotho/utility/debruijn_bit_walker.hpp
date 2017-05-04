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
#ifndef CLOTHO_DEBRUIJN_BIT_WALKER_HPP_
#define CLOTHO_DEBRUIJN_BIT_WALKER_HPP_

#include "clotho/utility/bit_masks.hpp"

namespace clotho {
namespace utility {

template < class T, unsigned char BlockSize >
struct debruijn_bit_walker_size;

template < class T >
struct debruijn_bit_walker_size< T, 8 > {
    static const unsigned int EMPTY = sizeof( T ) * 8;

    static unsigned int unset_next_index( T & b ) {
        if( b == 0 ) return EMPTY;
        
        unsigned int lsb = (unsigned int) b;
        lsb = LEAST_SIG_BIT( lsb );
        lsb = DEBRUIJNBIT_HASH_LOOKUP( lsb );

        b ^= bit_masks[ lsb ];

        return lsb;
    }
};

template < class T >
struct debruijn_bit_walker_size< T, 16 > {
    static const unsigned int EMPTY = sizeof( T ) * 8;
    static unsigned int unset_next_index( T & b ) {
        if( b == 0 ) return EMPTY;
        
        unsigned int lsb = (unsigned int) b;
        lsb = LEAST_SIG_BIT( lsb );
        lsb = DEBRUIJNBIT_HASH_LOOKUP( lsb );

        b ^= bit_masks[ lsb ];

        return lsb;
    }
};

template < class T >
struct debruijn_bit_walker_size< T, 32 > {
    static const unsigned int EMPTY = sizeof( T ) * 8;
    static unsigned int unset_next_index( T & b ) {
        if( b == 0 ) return EMPTY;
        
        unsigned int lsb = (unsigned int) b;
        lsb = LEAST_SIG_BIT( lsb );
        lsb = DEBRUIJNBIT_HASH_LOOKUP( lsb );

        b ^= bit_masks[ lsb ];

        return lsb;
    }

    static unsigned int next_and_shift( T & b ) {
        unsigned int lo = LEAST_SIG_BIT( b );
        unsigned int x = DEBRUIJNBIT_HASH_LOOKUP( lo );

        b >>= x;
        b ^= (T)1;

        return x;
    }
};

template < class T >
struct debruijn_bit_walker_size< T, 64 > {
    static const unsigned int EMPTY = sizeof( T ) * 8;
    // unset next lowest set bit
    //
    // return next set bit index
    // 
//    static unsigned int unset_next_index( T & b ) {
//        if( b == 0 ) return EMPTY;
//
//        unsigned int lsb = (unsigned int) b;
//
//        if( lsb == 0 ) {
//            lsb = (unsigned int)(b >> 32);
//            lsb = LEAST_SIG_BIT( lsb );
//            lsb = DEBRUIJNBIT_HASH_LOOKUP_HI( lsb );
//        } else {
//            lsb = LEAST_SIG_BIT( lsb );
//            lsb = DEBRUIJNBIT_HASH_LOOKUP( lsb );
//        }
//
////        b ^= bit_masks[ lsb ];
//        b &= ~(((T) 1) << lsb);
//
//        return lsb;
//    }

    static unsigned int next_and_shift( T & b ) {
        unsigned int offset = 0;
        unsigned int lo = (unsigned int) b;
        if( !lo ) {
            offset += 32;
            b >>= 32;

            lo = (unsigned int) b;
        }
        lo = LEAST_SIG_BIT( lo );
        unsigned int x = DEBRUIJNBIT_HASH_LOOKUP( lo );

        b >>= x;
        b ^= (T)1;

        return offset + x;
    }
};

template < class T >
struct debruijn_bit_walker : public debruijn_bit_walker_size< T, sizeof(T) * 8 > {
};

}   // namespace utility
}   // namespace clotho

#endif  // CLOTHO_DEBRUIJN_BIT_WALKER_HPP_
