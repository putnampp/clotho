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
#ifndef ALLELE_FREQUENCY_KERNELS_HPP_
#define ALLELE_FREQUENCY_KERNELS_HPP_

#include "clotho/cuda/data_spaces/sequence_space/device_sequence_space.hpp"
#include "clotho/cuda/data_spaces/basic_data_space.hpp"

template< class IntType, class CountType >
__global__ void count_alleles( device_sequence_space< IntType > * sequences, basic_data_space< CountType > * counts ) {
    unsigned int s_size = sequences->size;

    if( s_size == 0 ) return;
    typedef typename basic_data_space< CountType >::value_type   count_type;

    __shared__ count_type   s_data[ 1024 ];

    unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;
    unsigned int s_width = sequences->seq_width;

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int warp_id = (tid >> 5);
    unsigned int lane_id = (tid & 31);

    unsigned int warp_per_block = ((blockDim.x * blockDim.y) >> 5);
    unsigned int seq_offset = warp_per_block * s_width;

    typedef typename device_sequence_space< IntType >::int_type int_type;

    int_type * seqs = sequences->sequences;

    while( bid < s_width ) {    // true for all threads in a block
        unsigned int s_idx = warp_id * s_width + bid;
        count_type count = 0;

        while( s_idx < s_size ) {  // true for all threads in a warp
            int_type b = seqs[ s_idx ];

            count += ((b >> lane_id) & 1);
            
            s_idx += seq_offset;
        }
        __syncthreads();

        s_data[ tid ] = count;
        __syncthreads();
        // reduce threads
        //
        if( warp_id == 0 ) {
            unsigned int i = 32 + lane_id;
            while( i < (blockDim.x * blockDim.y) ) {
                count += s_data[ i ];
                i += 32;
            }

            counts->data[ bid * 32 + lane_id ] = count;
        }
        __syncthreads();

        bid += (gridDim.x * gridDim.y);
    }
}

/**
 * 
 * 1 thread per allele
 *
 */
template < class IntType >
__global__ void count_alleles( IntType * seqs, unsigned int * freq, unsigned int seq_count, unsigned int width, unsigned int allele_count ) {
    unsigned int all_idx = blockIdx.x * blockDim.x * blockDim.y + threadIdx.y * blockDim.x + threadIdx.x;

    unsigned int b_idx = blockIdx.x * blockDim.y + threadIdx.y;
    if( b_idx >= width ) {
        b_idx = width - 1;
    }
    __syncthreads();

    IntType mask = (1 << threadIdx.x);
    
    unsigned int N = 0;
    for( unsigned int i = 0; i < seq_count; ++i ) {
        IntType b = seqs[ i * width + b_idx ];

        if( b & mask ) {
            N += 1;
        }
    }
    __syncthreads();

    if( all_idx < allele_count ) {
        freq[ all_idx ] = N;
    }
}

#endif  // ALLELE_FREQUENCY_KERNELS_HPP_
