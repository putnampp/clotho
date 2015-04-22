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

// Popcount algorithm found on Hamming Weight wiki page:
// http://en.wikipedia.org/wiki/Hamming_weight
// (date 1/1/2015)
//
const uint64_t m1  = 0x5555555555555555; //binary: 0101...
const uint64_t m2  = 0x3333333333333333; //binary: 00110011..
const uint64_t m4  = 0x0f0f0f0f0f0f0f0f; //binary:  4 zeros,  4 ones ...
const uint64_t m8  = 0x00ff00ff00ff00ff; //binary:  8 zeros,  8 ones ...
const uint64_t m16 = 0x0000ffff0000ffff; //binary: 16 zeros, 16 ones ...
const uint64_t m32 = 0x00000000ffffffff; //binary: 32 zeros, 32 ones
const uint64_t hff = 0xffffffffffffffff; //binary: all ones
const uint64_t h01 = 0x0101010101010101; //the sum of 256 to the power of 0,1,2,3...

inline unsigned int popcount( unsigned int x ) {
    x -=((x >> 1) & (unsigned int) m1);
    x = (x & (unsigned int) m2) + ((x >> 2) & (unsigned int) m2);
    x = (x + (x >> 4)) & (unsigned int) m4;
    return (x * (unsigned int) h01) >> 24;
}

inline unsigned int popcount( unsigned long x ) {
    x -= (x >> 1) & m1;
    x = (x & m2) + ((x >> 2) & m2);
    x = (x + (x >> 4) ) & m4;
    return ( x * h01) >> 56;
}

#endif  // POPCOUNT_HPP_
