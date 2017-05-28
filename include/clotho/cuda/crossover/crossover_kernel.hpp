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
#ifndef CROSSOVER_KERNEL_HPP_
#define CROSSOVER_KERNEL_HPP_


/**
 *
 * 1 thread per sequence block
 * 1 block per 
 */
template < class IntType >
__global__ void crossover_kernel( IntType * parent, IntType * offspring, unsigned int * parent_indices, unsigned int parent_width, unsigned int offspring_width ) {
    assert( blockDim.x == 32 );

    unsigned int seq_block_idx = blockIdx.y * blockDim.x * blockDim.y + threadIdx.y * blockDim.x + threadIdx.x;

    unsigned int p_top = parent_indices[ blockIdx.x ];
    unsigned int p_bottom = p_top;

    if( p_bottom % 2 ) {
        p_bottom -= 1;
    } else {
        p_bottom += 1;
    }

    IntType top = 0, bottom = 0;

    if( seq_block_idx < parent_width ) {
        top = parents[ p_top * parent_width + seq_block_idx ];
        bottom = parents[ p_bottom * parent_width + seq_block_idx ];
    }
    __syncthreads();

    IntType mask = 0;

    if( seq_block_idx < offspring_width ) {
        mask = offspring[ blockIdx.x * offspring_width + seq_block_idx ];
    }
    __syncthreads();

    IntType off = (( top & ~mask ) | ( bottom & mask));

    if( seq_block_idx < offspring_width ) {
        offspring[ blockIdx.x * offspring_width + seq_block_idx ] = off;
    }
}

#endif  // CROSSOVER_KERNEL_HPP_
