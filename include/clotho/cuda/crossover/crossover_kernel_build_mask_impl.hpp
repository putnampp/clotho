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

    __shared__ int_type     event_count = 0;
    __shared__ real_type    rand_pool[ pool_type::MAX_EVENTS_PER_OFFSET ];

    // 1 sequence per block
    while( seqIdx < nSequences ) {
        int_type lo = events->offsets[ seqIdx ], hi = events->offsets[ seqIdx + 1 ];

        // true for all threads
        if( lo == hi ) {
            unsigned int id = tid;
            while( id < seq_width ) {
                buf[ seqIdx * seq_width + tid ] = 0;
                id += tpb;
            }
        } else {
            // move events from global memory into shared memory
            if( lo + tid < hi ) {
                rand_pool[tid] = events->event_pool[ lo + tid ];
            }
            __syncthreads();

            hi -= lo;

            unsigned int id = tid;
            unsigned int offset = seqIdx * seq_width + threadIdx.y;

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
                offset += wpb;
            }
        }

        seqIdx += bpg;
        __sycnthreads();
    }
}

#endif  // CROSSOVER_KERNEL_BUILD_MASK_IMPL_HPP_
