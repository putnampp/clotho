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
#ifndef CUDA_PAIRWISE_DIFFERENCE_KERNEL_HPP_
#define CUDA_PAIRWISE_DIFFERENCE_KERNEL_HPP_

#include "clotho/cuda/data_space/sequence_space/device_sequence_space.hpp"

#include "clotho/cuda/popcount_kernel.h"

struct pairwise_diff_stats {
    static const unsigned int BINS = 32;
    
    unsigned long long block_bin[ BINS ];

    unsigned long long count, total;
    double mean, stdev;
};

template < class IntType >
__global__ void pairwise_difference_kernel( device_sequence_space< IntType > * sequences, pairwise_diff_stats * stats ) {
    unsigned int R = sequences->seq_count;
    unsigned int C = sequences->seq_width;

    IntType * data = sequences->sequences;

    popcountGPU< IntType > pc;

    unsigned int bpg = (gridDim.x * gridDim.y); // blocks per grid

    unsigned long long count = R;
    count *= (count - 1);
    count <<= 1;

    unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    unsigned int s0 = bid * blockDim.y + threadIdx.y;

    unsigned int column_tail = (C % blockDim.x);
    unsigned int column_full_warps = (C / blockDim.x);
    column_full_warps *= blockDim.x;

    __shared__ unsigned int block_buffer[ 32 ]; // 32 == warps per block

    if(threadIdx.y == 0 ) {
        block_buffer[ threadIdx.x ] = 0;
    }
    __syncthreads();

    while( s0 < N ) {
        unsigned int s1 = s0 + 1;

        while( s1 < N ) {
            unsigned int s0_end = s0 * C + column_full_warps;
            unsigned int s0_start = s0 * C + threadIdx.x;

            unsigned int s1_offset = s1 * C + threadIdx.x;

            unsigned int tot = 0;
            while( s0_start < s0_end ) {    // all sequences same length so only need to iterate along first

                IntType s0_data = data[ s0_start ];
                IntType s1_data = data[ s1_start ];

                tot += pc.evalGPU( s0_data ^ s1_data );

                s0_start += blockDim.x;
                s1_start += blockDim.x;
            }
            __syncthreads();

            // handle tail
            if( column_tail ) { // true for all threads
                s0_end += column_tail;

                IntType s0_data = 0;
                IntType s1_data = 0;

                if( s0_start < s0_end ) {
                    s0_data = data[s0_start];
                    s1_data = data[s1_start];
                }
                __syncthreads();

                tot += pc.evalGPU( s0_data ^ s1_data );
            }
            __syncthreads();

            for( unsigned int i = 1; i < 32; i <<= 1 ) {
                unsigned int t = __shfl_up( tot, i );
                tot += ((threadIdx.x >= i) * t);
            }

            if( threadIdx.x == 31 ) {
                block_buffer[threadIdx.y] += tot;
            }
            __syncthreads();
        }
        s0 += bpg;
    }
}

#endif  // CUDA_PAIRWISE_DIFFERENCE_KERNEL_HPP_
