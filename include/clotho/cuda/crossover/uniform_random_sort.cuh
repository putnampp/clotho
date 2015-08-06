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
#ifndef UNIFORM_RANDOM_SORT_CUH_
#define UNIFORM_RANDOM_SORT_CUH_

#include <cuda.h>

template < unsigned int WarpSize >
struct uniform_random_sort {

    template < class RealType >
    __device__ void prefix_sum( RealType & r, unsigned int lane_id ) {
#pragma unroll
        for( unsigned int i = 1; i < WarpSize; i <<= 1 ) {
            RealType tmp = __shfl_up( r, i, WarpSize );
            if( lane_id >= i ) r += tmp;
        }
    }

    template < class RealType >
    __device__ RealType sort( volatile RealType * sData, RealType r, RealType pad, unsigned int tid ) {
        unsigned int lane_id = tid % WarpSize;
        unsigned int warp_id = tid / WarpSize;

        r = -log( r );

        prefix_sum( r, lane_id );

        if( lane_id == WarpSize - 1) {
            sData[ warp_id ] = r;
        }
        __syncthreads();

        RealType accum = sData[ lane_id ];
        prefix_sum( accum, lane_id );

        RealType tmp = __shfl( accum, warp_id - 1, WarpSize );
        if( warp_id > 0 ) { r += tmp; }

        tmp = __shfl( accum, WarpSize - 1, WarpSize );
        accum = tmp - log( pad );

        return r / accum;
    }
};

#endif  // UNIFORM_RANDOM_SORT_CUH_
