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
#ifndef CROSSOVER_KERNEL_BUILD_MASK_IMPL_HPP_
#define CROSSOVER_KERNEL_BUILD_MASK_IMPL_HPP_

#include "clotho/cuda/data_spaces/allele_space/device_allele_space.hpp"
#include "clotho/cuda/data_spaces/sequence_space/device_sequence_space.hpp"
#include "clotho/cuda/data_spaces/event_space/device_event_pool_space.hpp"

#include <cuda.h>

/**
 * This algorithm assumes that the event pool contains recombination events
 * have been generated elsewhere
 *
 */
template < class AlleleSpaceType, class RealType, class IntType >
__global__ void build_crossover_mask_kernel(  AlleleSpaceType * alleles,
    device_sequence_space< IntType > * buffer,
    device_event_pool_space< RealType, IntType > * events
    ) {

    typedef device_event_pool_space< RealType, IntType > pool_type;

    int_type bid = blockIdx.y * gridDim.x + blockIdx.x;
    int_type tid = threadIdx.y * blockDim.x + threadIdx.x;

    unsigned int tpb = (blockDim.x * blockDim.y);
    unsigned int bpg = (gridDim.x * gridDim.y);

    real_type * allele_list = alleles->locations;

    int_type * buf = buffer->buffer;
    int_type seq_width = buffer->seq_width;
    int_type nSequences = buffer->seq_count;

    unsigned int nAlleles = alleles->capacity;
    if( nAlleles == 0 ) return;

    __shared__ real_type    rand_pool[ pool_type::MAX_EVENTS_PER_OFFSET ];

    // 1 sequence per block
    while( bid < nSequences ) {
        int_type lo = events->offsets[ bid ], hi = events->offsets[ bid + 1 ];

        // true for all threads
        if( lo == hi ) {
            unsigned int id = bid * seq_width + tid;
            while( id < seq_width ) {
                buf[ id ] = 0;
                id += tpb;
            }
            __syncthreads();
        } else {
            // move events from global memory into shared memory
            if( lo + tid < hi ) {
                rand_pool[tid] = events->event_pool[ lo + tid ];
            }
            __syncthreads();

            hi -= lo;

            unsigned int id = tid;
            unsigned int offset = bid * seq_width + threadIdx.y;

            while( id < nAlleles ) {
                real_type pos = allele_list[ id ];

                unsigned int n = 0;
                for( unsigned int i = 0; i < hi; ++i ) {
                    n += ((rand_pool[ i ] < pos) ? 1 : 0);
                }
                __syncthreads();

                int_type mask = (1 << threadIdx.x);
                mask *= (n & 1);

                // scan right mask
                for( unsigned int i = 1; i < 32; i <<= 1 ) {
                    int_type _c = __shfl_up( mask, i );
                    mask |= ((threadIdx.x >= i ) *_c );
                }

                // last thread of each warp writes mask block to global memory
                if( threadIdx.x == 31 ) {
                    buf[ offset ] = mask;
                }
                __syncthreads();

                id += tpb;
                offset += blockDim.y;
            }
            __syncthreads();
        }

        bid += bpg;
    }
}

template < class RealType, class IntType >
__global__ void build_crossover_mask_kernel( RealType * locations, RealType * event_pool, unsigned int * event_dist, IntType * mask_buffer, unsigned int WIDTH, unsigned int ALLELE_COUNT ) {
    assert( blockDim.x == 32 );

    unsigned int seq_offset = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned int seq_idx = blockIdx.x * WIDTH + seq_offset;

    unsigned int event_start = event_dist[ blockIdx.x ];
    unsigned int event_end = event_dist[ blockIdx.x + 1 ];

    if( ALLELE_COUNT == 0 || event_start == event_end ) { // true for all threads
        if( threadIdx.x == 0 && seq_offset < WIDTH ) {
            mask_buffer[ seq_idx ] = 0;
        }
        __syncthreads();
    } else {
        unsigned int all_idx = blockIdx.y * blockDim.x * blockDim.y + threadIdx.y * blockDim.x + threadIdx.x;

        __shared__ RealType sEvents[ 32 ];
        if( threadIdx.y == 0 &&  event_start + threadIdx.x < event_end ) {
            sEvents[ threadIdx.x ] = event_pool[ threadIdx.x + event_start ]; 
        }
        __syncthreads();

        RealType loc = 1.0;

        if( all_idx < ALLELE_COUNT) {
            loc = locations[ all_idx ];
        }
        __syncthreads();

        IntType N = 0;
        while( event_start < event_end ) {
            RealType w = sEvents[ event_start++ ];
            N += (( loc < w ) ? 1 : 0 ); 
        }
        __syncthreads();

        IntType mask = (1 << threadIdx.x);
        mask *= (N & 1);

        for( unsigned int i = 1; i < 32; i <<= 1 ) {
            IntType _m = __shfl_up( mask, i );
            mask |= ((threadIdx.x >= i) * _m );
        }

        if( threadIdx.x == 31 && seq_offset < WIDTH ) {
            mask_buffer[ seq_idx ] = mask;
        }
        __syncthreads();
    } 
}

struct build_crossover_mask {
    template < class RealType, class IntType >
    static void execute( RealType * locations, RealType * event_pool, unsigned int * event_dist, IntType * mask_buffer, unsigned int seq_count, unsigned int seq_width, unsigned int allele_count ) {
        assert( locations != NULL );
        assert( event_pool != NULL );
        assert( event_dist != NULL );
        assert( mask_buffer != NULL );

        // 1 thread per device allele
        dim3 blocks( seq_count, 1, 1), threads( 1,1,1);

        unsigned int max_locations = seq_width * 32;
        assert( allele_count <= max_locations );

        if( max_locations > 1024 ) {
            blocks.y = max_locations / 1024;
            if( max_locations % 1024 ) {
                blocks.y += 1;
            }

            threads.x = 32;
            threads.y = 32;
        } else {
            threads.x = 32;
            threads.y = max_locations / 32;
            if( max_locations % 32 ) {
                threads.y += 1;
            }
            assert( threads.y <= 32);
        }

        std::cerr << "Build Mask - [ " << blocks.x << ", " << blocks.y << " ]; [ " << threads.x << ", " << threads.y << " ]" << std::endl;
        build_crossover_mask_kernel<<< blocks, threads >>>( locations, event_pool, event_dist, mask_buffer, seq_width, allele_count );
    }
};

#endif  // CROSSOVER_KERNEL_BUILD_MASK_IMPL_HPP_
