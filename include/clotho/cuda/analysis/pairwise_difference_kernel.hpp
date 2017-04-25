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

#include "clotho/cuda/analysis/pairwise_difference_stats.hpp"
#include "clotho/cuda/data_spaces/sequence_space/device_sequence_space.hpp"

#include "clotho/cuda/popcount_kernel.h"

__global__ void initialize_pairwise_diff_stats( pairwise_diff_stats * stats ) {
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int tpb = blockDim.x * blockDim.y;

    while( tid < pairwise_diff_stats::BINS ) {
        stats->block_bin[tid] = 0;
        tid += tpb;
    }
}

__global__ void finalize_pairwise_diff_stats( pairwise_diff_stats * stats ) {
    assert( blockDim.x * blockDim.y == 32 );

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    unsigned long long tot = stats->block_bin[tid];
    for( unsigned int i = 1; i < 32; i <<= 1 ) {
        unsigned long long t = __shfl_up( tot, i );
        tot += ((tid >= i ) ? t : 0);
    }

    if( tid == 31 ) {
        stats->total = tot;

        double d = (double)tot;
        double c = (double) stats->count;

        d /= c;

        stats->mean = d;
    }
}

template < class IntType >
__global__ void pairwise_difference_kernel( device_sequence_space< IntType > * sequences, basic_data_space< unsigned int > * sub_pop, pairwise_diff_stats * stats ) {

    unsigned int tpb = (blockDim.x * blockDim.y);

    assert( tpb % 32 == 0 );    // 32 == threads per warp; all warps are full

    unsigned int wpb = (tpb >> 5);  // warp per block

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int lane_id = (tid & 31);
    unsigned int warp_id = (tid >> 5);

    unsigned int bpg = (gridDim.x * gridDim.y); // blocks per grid
    unsigned int wpg = (wpb * bpg); // warps per grid

    unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;

    unsigned int R = sequences->seq_count;
    unsigned int C = sequences->seq_width;

    IntType * data = sequences->sequences;

    popcountGPU< IntType > pc;

    unsigned int N = sub_pop->size;
    unsigned int * sp = sub_pop->data;

    unsigned int column_tail = (C & 31);    // == % 32
//    unsigned int column_full_warps = (C & (~31));  // ==( / 32) * 32
    unsigned int column_full_warps = (C - column_tail);

    __shared__ unsigned long long block_buffer[ 32 ]; // 32 == max warps per block
    __shared__ unsigned int warp_index[32];

    unsigned int M = N - 1; // max index

    if( warp_id == 0 ) {
        block_buffer[ lane_id ] = 0;
        warp_index[ lane_id ] = M;
    }
    __syncthreads();

    // within each block verify that the sub-population is valid
    // NOTE consider moving this to a pre-processing step
    unsigned int s0 = tid;
    while( s0 < N ) {
        unsigned int idx = sp[s0];
        assert( idx < R );   
        s0 += tpb;
    }
    __syncthreads();


    s0 = 0;
    unsigned int s1 = bid * wpb + warp_id + 1;    // s1 = s0 + grid_warp_id + 1 =  warp_id + 1

    bool is_done = false;
    while( !is_done ) {

        while( s1 >= N  && s0 < M ) {
            ++s0;
            s1 -= (M - s0);
        }

        unsigned int s0_p = ((s0 >= M) ? M : s0);   // M is valid index in sub population list
        unsigned int s1_p = ((s1 >= M) ? M : s1);

        unsigned int s0_idx = sp[s0_p];
        unsigned int s1_idx = sp[s1_p];
        __syncthreads();

        unsigned int s0_end = s0_idx * C + column_full_warps;
        unsigned int s0_start = s0_idx * C + lane_id;

        unsigned int s1_start = s1_idx * C + lane_id;

        unsigned int tot = 0;
        while( s0_start < s0_end ) {    // all sequences same length so only need to iterate along first
            IntType s0_data = data[ s0_start ];
            IntType s1_data = data[ s1_start ];

            tot += pc.evalGPU( s0_data ^ s1_data );

            s0_start += 32;
            s1_start += 32;
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
            tot += ((lane_id >= i) * t);
        }

        if( lane_id == 31 ) {
            block_buffer[warp_id] += tot;
            warp_index[warp_id] =  s0;
        }
        __syncthreads();

        // find the minimum s0 within the block
        unsigned int tmp = warp_index[lane_id];
        for( unsigned int i = 1; i < 32; i <<= 1 ) {
            unsigned int t = __shfl_up(tmp, i);
            tmp = ((tmp < t) ? tmp : t);
        }

        unsigned int t = __shfl(tmp, 31);

        is_done = (t >= M);   // done when all warps have reach s0 >= M
        s1 += wpg;
    }

    unsigned long long tot = block_buffer[ lane_id ];
    for( unsigned int i = 1; i < 32; i <<= 1 ) {
        unsigned long long t = __shfl_up( tot, i );
        tot += (( lane_id >= i ) ? t : 0);
    }

    if( warp_id == 0 && lane_id == 31 ) {
        stats->block_bin[ bid ] = tot;

        if(bid == 0) {
            unsigned long long count = N;
            count *= (count - 1);
            count >>= 1;
            stats->count = count;
        }
    }
}

template < class IntType >
__global__ void pairwise_difference_kernel2( device_sequence_space< IntType > * sequences, basic_data_space< unsigned int > * sub_pop, pairwise_diff_stats * stats ) {

    unsigned int tpb = (blockDim.x * blockDim.y);

    assert( tpb % 32 == 0 );    // 32 == threads per warp; all warps are full

    //unsigned int wpb = (tpb >> 5);  // warp per block

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int lane_id = (tid & 31);
    unsigned int warp_id = (tid >> 5);

    unsigned int bpg = (gridDim.x * gridDim.y); // blocks per grid
//    unsigned int wpg = (wpb * bpg); // warps per grid

    unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;

    unsigned int R = sequences->seq_count;
    unsigned int C = sequences->seq_width;

    IntType * data = sequences->sequences;

    popcountGPU< IntType > pc;

    unsigned int N = sub_pop->size;
    unsigned int * sp = sub_pop->data;

    unsigned int column_tail = (C % tpb);
    unsigned int column_full_warps = (C - column_tail);

    __shared__ unsigned long long block_buffer[ 32 ]; // 32 == max warps per block

    unsigned int M = N - 1; // max index

    if(bid == 0 && tid == 0) {
        unsigned long long count = N;
        count *= (count - 1);
        count >>= 1;
        stats->count = count;   // count == (N * (N -1)) / 2
    }
    __syncthreads();

    if( warp_id == 0 ) {
        block_buffer[ lane_id ] = 0;
    }
    __syncthreads();

    // within each block verify that the sub-population is valid
    // NOTE consider moving this to a pre-processing step
    unsigned int s0 = tid;
    while( s0 < N ) {
        unsigned int idx = sp[s0];
        assert( idx < R );   
        s0 += tpb;
    }
    __syncthreads();


    s0 = 0;
    unsigned int s1 = bid + 1;

    while( s0 < M ) {
        while( s1 >= N  && s0 < M ) {
            ++s0;
            s1 -= (M - s0);
        }

        unsigned int s0_p = ((s0 >= M) ? M : s0);   // M is valid index in sub population list
        unsigned int s1_p = ((s1 >= M) ? M : s1);

        unsigned int s0_idx = sp[s0_p]; // get sub-population index corresponding to each index
        unsigned int s1_idx = sp[s1_p];
        __syncthreads();

        unsigned int s0_end = s0_idx * C + column_full_warps;
        unsigned int s0_start = s0_idx * C + tid;

        unsigned int s1_start = s1_idx * C + tid;

        unsigned int tot = 0;
        while( s0_start < s0_end ) {    // all sequences same length so only need to iterate along first
            IntType s0_data = data[ s0_start ];
            IntType s1_data = data[ s1_start ];

            tot += pc.evalGPU( s0_data ^ s1_data );

            s0_start += tpb;
            s1_start += tpb;
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

        // scan right of lane totals
        for( unsigned int i = 1; i < 32; i <<= 1 ) {
            unsigned int t = __shfl_up( tot, i );
            tot += ((lane_id >= i) ? t : 0 );
        }

        if( lane_id == 31 ) {
            block_buffer[warp_id] += tot;
        }
        __syncthreads();

        s1 += bpg;
    }

    if( warp_id == 0 ) {
        unsigned long long tot = block_buffer[ lane_id ];
        for( unsigned int i = 1; i < 32; i <<= 1 ) {
            unsigned long long t = __shfl_up( tot, i );
            tot += (( lane_id >= i ) ? t : 0);
        }

        if( lane_id == 31 ) {
            stats->block_bin[ bid ] = tot;
        }
    }
}

#endif  // CUDA_PAIRWISE_DIFFERENCE_KERNEL_HPP_
