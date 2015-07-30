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
#ifndef WARP_SORT_32_HPP_
#define WARP_SORT_32_HPP_

#include "clotho/cuda/warp_sort/warp_sort_def.hpp"

namespace clotho {
namespace cuda {

template <>
struct warp_sort< 32 > {
    static const unsigned int MAX_BITONIC_ROUNDS = 5;

    template < class K >
    __device__ void sort( K & k, unsigned int round ) {
        unsigned int round_mask = ((round == MAX_BITONIC_ROUNDS) ? 0 : ((threadIdx.x >> round) & 0x01));
        unsigned int xor_mask = (1 << (round - 1));

        while( round-- ) {
            unsigned int dir = ((round_mask ^ (threadIdx.x >> round)) & 0x01);
            K sw_k = __shfl_xor( k, xor_mask );
            if( k < sw_k == dir ) {
                k = sw_k;
            }

            xor_mask >>= 1;
        }
    }

    template < class K, class V >
    __device__ void sort( K & k, V & v, unsigned int round ) {
        unsigned int round_mask = ((round == MAX_BITONIC_ROUNDS) ? 0 : ((threadIdx.x >> round) & 0x01));
        unsigned int xor_mask = (1 << (round - 1));

        while( round-- ) {
            unsigned int dir = ((round_mask ^ (threadIdx.x >> round)) & 0x01);
            K sw_k = __shfl_xor( k, xor_mask );
            V sw_v = __shfl_xor( v, xor_mask );
            if( k < sw_k == dir ) {
                k = sw_k;
                v = sw_v;
            }

            xor_mask >>= 1;
        }
    }
};

}   // namespace cuda {
}   // namespace clotho {

#endif  // WARP_SORT_32_HPP_
