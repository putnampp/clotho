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
#ifndef CLOTHO_CUDA_DISCRETE_DISTRIBUTION_HPP_
#define CLOTHO_CUDA_DISCRETE_DISTRIBUTION_HPP_

#include <cuda.h>

template < class RealType >
__global__ void normalize_kernel( RealType * in, RealType * out, unsigned int N ) {
    typedef RealType real_type;

    // set block relative thread id
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int warp_id = (tid >> 5);  // assumes 32 threads/warp
    unsigned int lane_id = (tid & 31);  // assumes 32 threads/warp

    unsigned int tcount = blockDim.x * blockDim.y;

    real_type _s = 0.;

    while( tid < N ) {
        real_type v = in[tid];
        _s += v;
        tid += tcount;
    }
    __syncthreads();

    // reset thread id
    tid = threadIdx.y * blockDim.x + threadIdx.x;

    for( unsigned int i = 1; i < 32; i <<= 1 ) {
        real_type s = __shfl_up( _s, i );
        _s += ((real_type)( lane_id >= i )) * s;
    }

    if( tcount > 32 ) { // true/false for all threads in block
        __shared__ real_type buffer[ 32 ];  // max of 32 warps/block

        if( warp_id == 0 ) {    // use the threads in warp_0 to clear buffer 
            buffer[lane_id] = 0.;
        }
        __syncthreads();

        buffer[ warp_id ] = _s;
        __syncthreads();

        _s = buffer[lane_id];

        for( unsigned int i = 1; i < 32; i <<= 1 ) {
            real_type s = __shfl_up( _s, i );
            _s += ((real_type)( lane_id >= i )) * s;
        }
    }
 
    // every warp computes the total summation into thread 31 of warp   
    real_type sum = __shfl( _s, 31 );

    while( tid < N ) {
        real_type v= in[tid];
        v /= sum;

        out[tid] = v;
        tid += tcount;
    }
}

#endif  // CLOTHO_CUDA_DISCRETE_DISTRIBUTION_HPP_
