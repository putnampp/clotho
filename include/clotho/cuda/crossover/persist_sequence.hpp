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
#ifndef PERSIST_SEQUENCE_HPP_
#define PERSIST_SEQUENCE_HPP_

#include <cuda.h>

template < class IntType >
__device__ void persist_mask( IntType mask, IntType warp_id, IntType lane_id, IntType * seq ) {
    // collapse masks to single crossover mask per warp
    // mask will exist in lane 0 for all warps
#pragma unroll
    for( unsigned int j = 1; j < 32; j <<= 1 ) {
        unsigned int tmp = __shfl_down( mask, j );
        mask |= ((!( lane_id & (( j << 1) - 1))) * tmp);
    }

    // use single thread per warp to write/store
    // crossover mask to global memory
    if( lane_id == 0) {
        seq[ warp_id ] = mask;
    }
}

template < class IntType >
__device__ void persist_mask_unrolled( IntType mask, IntType warp_id, IntType lane_id, IntType * seq ) {
    // collapse masks to single crossover mask per warp
    // mask will exist in lane 0 for all warps

    unsigned int tmp = __shfl_down( mask, 1 );
    mask |= ((!( lane_id & 1)) * tmp);
//    mask |= ((( lane_id & 1) == 0) ? tmp : 0 );

    tmp = __shfl_down( mask, 2);
    mask |= ((!( lane_id & 3)) * tmp);
//    mask |= ((( lane_id & 3) == 0) ? tmp : 0 );
    
    tmp = __shfl_down( mask, 4);
    mask |= ((!( lane_id & 7)) * tmp);
//    mask |= ((( lane_id & 7) == 0) ? tmp : 0 );

    tmp = __shfl_down( mask, 8);
    mask |= ((!( lane_id & 15)) * tmp);
//    mask |= ((( lane_id & 15) == 0) ? tmp : 0 );

    tmp = __shfl_down( mask, 16);
    mask |= ((!( lane_id & 31)) * tmp);
//    mask |= ((( lane_id & 31) == 0) ? tmp : 0 );

    // use single thread per warp to write/store
    // crossover mask to global memory
    if( lane_id == 0) {
        seq[ warp_id ] = mask;
    }
}

template < class IntType >
__device__ void persist_mask( unsigned char mask, IntType warp_id, IntType lane_id, unsigned char * seq ) {
    seq[ warp_id * 32 + lane_id ] = mask;
}
#endif  // PERSIST_SEQUENCE_HPP_
