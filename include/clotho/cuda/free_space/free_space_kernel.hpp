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
#ifndef FREE_SPACE_KERNEL_HPP_
#define FREE_SPACE_KERNEL_HPP_

#include "clotho/cuda/free_space/offset_enum.hpp"

template < class IntType >
__global__ void evaluate_free_space( IntType * pop, IntType * fspace, unsigned int seq_count, unsigned int width ) {

    unsigned int offset = blockIdx.x * blockDim.x * blockDim.y + threadIdx.y * blockDim.x + threadIdx.x;

    if( offset >= width ) {
        offset = width - 1;
    }
    __syncthreads();

    IntType var = 0;
    IntType fixed = ~var;

    for( unsigned int i = 0; i < seq_count; ++i ) {
        IntType b = pop[ i * width + offset ];

        var |= b;
        fixed &= b;
    }

    __syncthreads();

    offset = blockIdx.x * blockDim.x * blockDim.y + threadIdx.y * blockDim.x + threadIdx.x;
    if( offset < width ) {
        fspace[ offset + LOST_OFFSET * width ] = ~var;
        fspace[ offset + FIXED_OFFSET *width ] = fixed;
        fspace[ offset + FREE_OFFSET * width ] = ~(var | fixed);
    }

    __syncthreads();
}

template < class IntType >
__global__ void remove_fixed( IntType * pop, IntType * fspace, unsigned int seq_count, unsigned int width ) { 

    unsigned int offset = blockIdx.x * blockDim.x * blockDim.y + threadIdx.y * blockDim.x + threadIdx.x;

    IntType fixed = 0;
    if( offset < width ) {
        fixed = fspace[ offset + FIXED_OFFSET * width ];
    }
    __syncthreads();

    for( unsigned int i = 0; i < seq_count; ++i ) {
        if( offset < width ) {
            IntType b = pop[ i * width + offset ];

            assert( (b & fixed) == fixed);
            b &= ~fixed;

            pop[ i * width + offset ] = b;
        }
        __syncthreads();
    }
}

#endif  // FREE_SPACE_KERNEL_HPP_
