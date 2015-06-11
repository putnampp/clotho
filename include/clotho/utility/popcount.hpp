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
#ifndef POPCOUNT_HPP_
#define POPCOUNT_HPP_

/*
// Popcount algorithm found on Hamming Weight wiki page:
// http://en.wikipedia.org/wiki/Hamming_weight
// (date 1/1/2015)
//
const unsigned long long M1  = 0x5555555555555555; //binary: 0101...
const unsigned long long M2  = 0x3333333333333333; //binary: 00110011..
const unsigned long long M4  = 0x0f0f0f0f0f0f0f0f; //binary:  4 zeros,  4 ones ...
const unsigned long long M8  = 0x00ff00ff00ff00ff; //binary:  8 zeros,  8 ones ...
const unsigned long long M16 = 0x0000ffff0000ffff; //binary: 16 zeros, 16 ones ...
const unsigned long long M32 = 0x00000000ffffffff; //binary: 32 zeros, 32 ones
const unsigned long long HFF = 0xffffffffffffffff; //binary: all ones
const unsigned long long H01 = 0x0101010101010101; //the sum of 256 to the power of 0,1,2,3...

inline unsigned int popcount( unsigned int x ) {
    x -=((x >> 1) & (unsigned int) M1);
    x = (x & (unsigned int) M2) + ((x >> 2) & (unsigned int) M2);
    x = (x + (x >> 4)) & (unsigned int) M4;
    return (x * (unsigned int) H01) >> 24;
}

inline unsigned long long popcount( unsigned long long x ) {
    x -= (x >> 1) & M1;
    x = (x & M2) + ((x >> 2) & M2);
    x = (x + (x >> 4) ) & M4;
    return ( x * H01) >> 56;
}*/

template < unsigned int WordWidth >
struct popcount_constants;

/**
 *  Byte Lookup Table Based popcount
 *  Allocates 1K (256 * 4 bytes) of memory
 *
 *  (6/11/15)
 *  Clever macro construction found at:
 *  http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetTable
 *
 */
template < >
struct popcount_constants< 1 > {
    static const unsigned int lookup[ 256 ]; 

    template < class T >
    static unsigned int eval( T x ) {
        return lookup[ (unsigned char) x ];
    }
};

template < >
struct popcount_constants< 2 > {
    template < class T >
    static unsigned int eval( T t ) {
        return eval( (unsigned short) t );
    }

    static unsigned int eval( unsigned short x ) {
        typedef popcount_constants< 1 > base_type;
        return (base_type::lookup[ x & 0x00FF ] + base_type::lookup[x >> 8]);
    }
};

template < >
struct popcount_constants< 4 > {
    static const unsigned int M1  = 0x55555555; //binary: 0101...
    static const unsigned int M2  = 0x33333333; //binary: 00110011..
    static const unsigned int M4  = 0x0f0f0f0f; //binary:  4 zeros,  4 ones ...
    static const unsigned int M8  = 0x00ff00ff; //binary:  8 zeros,  8 ones ...
    static const unsigned int M16 = 0x0000ffff; //binary: 16 zeros, 16 ones ...
    static const unsigned int HFF = 0xffffffff; //binary: all ones
    static const unsigned int H01 = 0x01010101; //the sum of 256 to the power of 0,1,2,3...
    static const unsigned long long DSHIFT = 24;

    template < class T >
    static unsigned int eval( T x ) {
        return eval( (unsigned int) x );
    }

    static unsigned int eval( unsigned int x ) {
        x -=((x >> 1) & M1);
        x = (x & M2) + ((x >> 2) & M2);
        x = (x + (x >> 4)) & M4;
        return (x * H01) >> DSHIFT;
    }
};

template < >
struct popcount_constants< 8 > {
    static const unsigned long long M1  = 0x5555555555555555; //binary: 0101...
    static const unsigned long long M2  = 0x3333333333333333; //binary: 00110011..
    static const unsigned long long M4  = 0x0f0f0f0f0f0f0f0f; //binary:  4 zeros,  4 ones ...
    static const unsigned long long M8  = 0x00ff00ff00ff00ff; //binary:  8 zeros,  8 ones ...
    static const unsigned long long M16 = 0x0000ffff0000ffff; //binary: 16 zeros, 16 ones ...
    static const unsigned long long M32 = 0x00000000ffffffff; //binary: 32 zeros, 32 ones
    static const unsigned long long HFF = 0xffffffffffffffff; //binary: all ones
    static const unsigned long long H01 = 0x0101010101010101; //the sum of 256 to the power of 0,1,2,3...
    static const unsigned long long DSHIFT = 56;

    template < class T >
    static unsigned int eval( T t ) {
        return eval( (unsigned long long) t );
    }

    static unsigned int eval( unsigned long long x ) {
        x -= (x >> 1) & M1;
        x = (x & M2) + ((x >> 2) & M2);
        x = (x + (x >> 4) ) & M4;
        return ( x * H01) >> DSHIFT;
    }
};

template < class T >
unsigned int popcount( T x ) {
    typedef popcount_constants< sizeof(T) > base_type;
    return base_type::eval( x );
}

#endif  // POPCOUNT_HPP_
