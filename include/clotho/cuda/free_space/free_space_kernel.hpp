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

/**
 * 
 * 1 thread per sequence block
 * 
 */
template < class IntType >
__global__ void evaluate_free_space_kernel( IntType * pop, IntType * fspace, unsigned int seq_count, unsigned int width ) {

    // blockDim.x * blockDim.y == sequence blocks per batch
    // threadIdx.y * blockDim.x + threadIdx.x == sequence block index

    unsigned int block_idx = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int seq_offset = blockIdx.x * blockDim.x * blockDim.y + block_idx;

    IntType var = 0;
    IntType fixed = (~var);

    if( seq_offset >= width ) {
        seq_offset = width - 1;
        fixed = 0;  // no fixed 
    }
    __syncthreads();

    for( unsigned int i = 0; i < seq_count; ++i ) {
        IntType b = pop[ i * width + seq_offset ];

        var = (var | b);
        fixed = (fixed & b);
    }
    __syncthreads();

    seq_offset = blockIdx.x * blockDim.x * blockDim.y + block_idx;
    if( seq_offset < width ) {
        fspace[ seq_offset + 0 * width ] = fixed;
        fspace[ seq_offset + 1 * width ] = (~var);
        fspace[ seq_offset + 2 * width ] = (~(var | fixed));
    }

    __syncthreads();
}

struct evaluate_free_space {
    template < class IntType >
    static void execute( IntType * pop, IntType * fspace, unsigned int seq_count, unsigned int width ){
        assert( pop != NULL );
        assert( fspace != NULL );
        assert( width % 32 == 0 );

        dim3 blocks( 1,1,1), threads( 1,1,1 );
        if( width > 1024 ) {
            blocks.x = (width / 1024) + ((width % 1024) ? 1 : 0);
            threads.x = 32;
            threads.y = 32;
        } else {
            threads.x = 32;
            threads.y = width / 32;
        }
        std::cerr << "Free Space - [ " << blocks.x << ", " << blocks.y << " ]; [ " << threads.x << ", " << threads.y << " ]" << std::endl;
        evaluate_free_space_kernel<<< blocks, threads >>>( pop, fspace, seq_count, width );
    }
};


/**
 *
 * 1 thread per sequence block
 *
 */
template < class IntType >
__global__ void remove_fixed_kernel( IntType * pop, IntType * fspace, unsigned int width ) { 

    unsigned int offset = blockIdx.y * blockDim.x * blockDim.y + threadIdx.y * blockDim.x + threadIdx.x;

    if( offset < width ) {
        IntType fixed = fspace[ offset + FIXED_OFFSET * width ];
        
        IntType b = pop[ blockIdx.x * width + offset ];

        assert( (b & fixed) == fixed);
        b = (b & (~fixed));

        pop[ blockIdx.x * width + offset ] = b;
    }
}

struct remove_fixed {

    template < class IntType > 
    static void execute( IntType * pop, IntType * fspace, unsigned int seq_count, unsigned int seq_width ) {
        assert( pop != NULL );
        assert( fspace != NULL );
        assert( seq_width % 32 == 0 );

        dim3 blocks( seq_count,1,1), threads( 1,1,1 );
        if( seq_width > 1024 ) {
            blocks.y = (seq_width / 1024) + ((seq_width % 1024) ? 1 : 0);
            threads.x = 32;
            threads.y = 32;
        } else {
            threads.x = 32;
            threads.y = seq_width / 32;
        }
        std::cerr << "Remove Fixed - [ " << blocks.x << ", " << blocks.y << " ]; [ " << threads.x << ", " << threads.y << " ]" << std::endl;
        remove_fixed_kernel<<< blocks, threads >>>( pop, fspace, seq_width );

    }
};

#endif  // FREE_SPACE_KERNEL_HPP_
