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
#ifndef CROSSOVER_MATRIX_5_CUH_
#define CROSSOVER_MATRIX_5_CUH_

#include "clotho/cuda/crossover/uniform_random_sort.cuh"
#include "clotho/cuda/crossover/crossover_matrix_def.hpp"

#include <cuda.h>

template <>
const unsigned int crossover< 5 >::ALLELE_PER_INT = 32;

template < >
const unsigned int crossover< 5 >::MAX_EVENTS = ((crossover< 5 >::comp_cap_type::MAX_CONSTANT_MEMORY / sizeof( crossover< 5 >::event_count_type )) >> 1);

template < class RealType, class CC >
__global__ void order_random( RealType * in_list, RealType * out_list, size_t N ) {
    if( N > 1024 ) return;

    typedef RealType real_type;

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    uniform_random_sort< CC::WARP_SIZE > order;
    __shared__ real_type sBuffer[ 1025 ];

    if( tid <= N ) {
        sBuffer[ tid ] = in_list[ tid ];
    } else {
        sBuffer[ tid ] = 1.0;
    }
    __syncthreads();

    real_type pad = sBuffer[N];
    real_type r = sBuffer[tid];
    if( tid == N ) {
        r = 1.0;
    }
    __syncthreads();

    real_type rpt = order.sort( sBuffer, r, pad, tid );
    if( tid < N ) {
        out_list[ tid ] = rpt;
    }
}

template < class RealType >
__global__ void select_randoms( RealType * rand_pool, unsigned int offset, unsigned int N, RealType * o_list ) {
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    if( tid < N ) {
        o_list[ tid ] = rand_pool[ offset + tid ];
    }
}

#ifdef USE_CONSTANT_EVENT_LIST
#define EVT_LIST g_cEvents

__constant__ crossover< 5 >::event_count_type EVT_LIST[ crossover< 5 >::MAX_EVENTS + 1 ];

template < class RealType, class CC >
__global__ void crossover_matrix_5( RealType * rand_pool
                                  , RealType * alleles
                                  , unsigned int * sequences
                                  , unsigned int sequence_width
                                  , unsigned int round_offset ) {
#else
#define EVT_LIST event_list

template < class RealType, class CC >
__global__ void crossover_matrix_5( RealType * rand_pool
                                  , RealType * alleles
                                  , unsigned int * EVT_LIST
                                  , unsigned int * sequences
                                  , unsigned int sequence_width
                                  , unsigned int round_offset ) {
#endif  // USE_CONSTANT_EVENT_LIST

    typedef RealType AlleleType;
    typedef unsigned int CountType;

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    CountType eStart = EVT_LIST[blockIdx.y], eEnd = EVT_LIST[blockIdx.y + 1];

    // if there are no recombination events for this thread block/sequence
    if( eStart >= eEnd ) {  // will be true or false for all threads in the block
        if( threadIdx.y == 0 ) {
            sequences[ blockIdx.y * sequence_width + blockIdx.x * blockDim.x + threadIdx.x ] = 0;
        }
        return;
    }

    // casting to int to allow for negative indices
    int min = 0, max = (eEnd - eStart) - 1;

    if( max >= CC::THREADS_PER_BLOCK ) {
        return;
    }

    // load an allele into a register
    AlleleType allele = alleles[ blockIdx.x * CC::THREADS_PER_BLOCK + tid ];
    __syncthreads();

    uniform_random_sort< CC::WARP_SIZE > order;

    __shared__ RealType recomb_points[ CC::THREADS_PER_BLOCK + 1 ];

    // grab an extra random number for padding of final element
    if( eStart + tid <= eEnd ) {
        recomb_points[ tid ] = rand_pool[ round_offset + blockIdx.y + eStart + tid ];
    } else {
        recomb_points[ tid ] = 1.0;
    }
    __syncthreads();

    RealType rpt = ((eStart + tid == eEnd ) ? 1.0 : recomb_points[tid]);
    RealType pad = recomb_points[ eEnd - eStart ];
    __syncthreads();
    
    rpt = order.sort( recomb_points, rpt, pad, tid );
    __syncthreads();

    recomb_points[ tid ] = rpt;
    __syncthreads();

    // binary search of recombination point list
    while( min <= max ) {
        int mid = ((max - min) / 2) + min;
        AlleleType tmp = recomb_points[ mid ];

        if( tmp < allele ) {
            min = mid + 1;
        } else if( tmp > allele ) {
            max = mid - 1;
        } else {
            // allele occurs at a recombination point
            min = mid;
            break;
        }
    }
    __syncthreads();
    // min contains count of preceeding recombination events

    unsigned int mask = ((min & 1) * (1 << threadIdx.x));

#pragma unroll
    for( unsigned int i = 1; i < CC::WARP_SIZE; i <<= 1) {
        unsigned int m = __shfl_down(mask, i, CC::WARP_SIZE);
        mask |= ( (!(threadIdx.x & ((i << 1) - 1))) * m);
    }

    if(threadIdx.x == 0 ) {
        sequences[ blockIdx.y * sequence_width + blockIdx.x * blockDim.x + threadIdx.y ] = mask;
    }
}

template <>
void crossover< 5 >::operator()(  real_type         * rand_pool
                            , allele_type       * allele_list
                            , event_count_type  * evt_list
                            , int_type          * sequences
                            , size_type nSequences
                            , size_type nAlleles
                            , size_type sequence_width  ) {
    assert( nAlleles % comp_cap_type::THREADS_PER_BLOCK == 0 );
    assert( sequence_width >= nAlleles / comp_cap_type::WARP_SIZE );

    unsigned int round_offset = 0;
    while( round_offset < nSequences ) {
        unsigned int N = (nSequences - round_offset);
        if( N > MAX_EVENTS ) { N = MAX_EVENTS; }

        // execution configuration:
        // blockIdx.x = group of THREAD_PER_BLOCK alleles (1 column tile per 1024 alleles)
        // blockIdx.y = offset of offspring sequence (1 - row tile per offspring sequence)
        // blockIdx.z = 1; unused
        //
        // 1 allele/thread therefore
        // threadIdx.x = warp offset
        // threadIdx.y = warps per block
        // threadIdx.z = 1; unused
        dim3 blocks( nAlleles / comp_cap_type::THREADS_PER_BLOCK, N, 1)
            , threads( comp_cap_type::WARP_SIZE, comp_cap_type::WARP_PER_BLOCK, 1);

#ifdef USE_CONSTANT_EVENT_LIST
        cudaMemcpyToSymbol( EVT_LIST, evt_list, (N + 1) * sizeof( event_count_type ), 0, cudaMemcpyDeviceToDevice );

        // assumes the execution configuration specifies sequences to generate
        // allele tile offset determined from execution configuration
        crossover_matrix_5< real_type, comp_cap_type ><<< blocks, threads >>>( rand_pool, allele_list, sequences, sequence_width, round_offset );
#else
        crossover_matrix_5< real_type, comp_cap_type ><<< blocks, threads >>>( rand_pool, allele_list, evt_list, sequences, sequence_width, round_offset );
#endif  // USE_CONSTANT_EVENT_LIST

        cudaDeviceSynchronize();

        evt_list += N;
        round_offset += N;
        sequences += (N * sequence_width);
    }
}

#endif  // CROSSOVER_MATRIX_5_CUH_
