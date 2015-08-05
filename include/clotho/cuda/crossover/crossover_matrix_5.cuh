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

#include "clotho/cuda/compute_capability.hpp"

#include <cuda.h>
#include <iostream>

class crossover {
public:
    typedef double                       real_type;
    typedef double                       allele_type;
    typedef unsigned int                event_count_type;
    typedef unsigned int                int_type;
    typedef unsigned int                size_type;
    typedef compute_capability< 3, 0 >  comp_cap_type;

    static const unsigned int MAX_EVENTS = ((comp_cap_type::MAX_CONSTANT_MEMORY / sizeof( event_count_type )) >> 1);    // use half of the constant space for event counts

    void operator()(  real_type         * rand_pool
                    , allele_type       * allele_list
                    , event_count_type  * event_list
                    , int_type          * sequences
                    , size_type nSequences
                    , size_type nAlleles
                    , size_type sequence_width );
 
    virtual ~crossover() {}
};


template < class RealType, class CC >
__global__ void order_random( RealType * in_list, RealType * out_list, size_t N ) {
    if( N > 1024 ) return;

    typedef RealType real_type;

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    __shared__ real_type sBuffer[ 1024 ];
    real_type rpt = 1.0;

    if( tid < N ) {
        rpt = in_list[ tid ];
    }
    __syncthreads();

    // order the list of recombination points
    rpt = -log( rpt );

    // scan
    for( unsigned int i = 1; i < CC::WARP_SIZE; i <<= 1 ) {
        real_type tmp = __shfl_up( rpt, i, CC::WARP_SIZE );
        if( threadIdx.x >= i ) rpt += tmp;
    }

    // share partial sums with other warps in block
    if( threadIdx.x == CC::WARP_SIZE - 1 ) {
        sBuffer[ threadIdx.y ] = rpt;
    }
    __syncthreads();

    real_type accum = sBuffer[threadIdx.x];
    for( unsigned int i = 1; i < CC::WARP_SIZE; i <<= 1 ) {
        real_type tmp = __shfl_up( accum, i, CC::WARP_SIZE );
        if( threadIdx.x >= i ) accum += tmp;
    }
    
    if( threadIdx.y > 0 ) {
        real_type tmp = __shfl( accum, threadIdx.y - 1, CC::WARP_SIZE );
        rpt += tmp;
    }
    __syncthreads();

    real_type tmp = __shfl( accum, CC::WARP_SIZE - 1, CC::WARP_SIZE );
    accum = tmp - log(0.1); // technically should be log of an additional random value (not a constant)

    rpt /= accum;

    if( tid < N ) {
        out_list[tid] = rpt;
    }
}

__constant__ crossover::event_count_type g_cEvents[ crossover::MAX_EVENTS + 1 ];

template < class RealType, class CC >
__global__ void crossover_matrix( RealType * rand_pool
                                  , RealType * alleles
                                  , unsigned int * sequences
                                  , unsigned int sequence_width ) {

    typedef RealType AlleleType;
    typedef unsigned int CountType;

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    CountType eStart = g_cEvents[blockIdx.y], eEnd = g_cEvents[blockIdx.y + 1];

    // if there are no recombination events for this thread block/sequence
    if( eStart >= eEnd ) {  // will be true or false for all threads in the block
        if( threadIdx.y == 0 ) {
            sequences[ blockIdx.y * sequence_width + blockIdx.x * blockDim.x + threadIdx.x ] = 0;
        }
        return;
    }

    __shared__ RealType recomb_points[ CC::THREADS_PER_BLOCK ];
    RealType rpt = 1.0;

    if( eStart + tid < eEnd ) {
        rpt = rand_pool[ eStart + tid ];
    }
    __syncthreads();

    // order the list of recombination points
    rpt = -log( rpt );

    // scan
    for( unsigned int i = 1; i < CC::WARP_SIZE; i <<= 1 ) {
        RealType tmp = __shfl_up( rpt, i, CC::WARP_SIZE );
        if( threadIdx.x >= i ) rpt += tmp;
    }

    // share partial sums with other warps in block
    if( threadIdx.x == CC::WARP_SIZE - 1 ) {
        recomb_points[ threadIdx.y ] = rpt;
    }
    __syncthreads();

    RealType accum = recomb_points[threadIdx.x];
    for( unsigned int i = 1; i < CC::WARP_SIZE; i <<= 1 ) {
        RealType tmp = __shfl_up( accum, i, CC::WARP_SIZE );
        if( threadIdx.x >= i ) accum += tmp;
    }
    
    if( threadIdx.y > 0 ) {
        RealType tmp = __shfl( accum, threadIdx.y - 1, CC::WARP_SIZE );
        rpt += tmp;
    }
    __syncthreads();

    RealType tmp = __shfl( accum, CC::WARP_SIZE - 1, CC::WARP_SIZE );
    accum = tmp - log(0.1); // technically should be log of an additional random value (not a constant)

    rpt /= accum;
    recomb_points[tid] = rpt;
    // at this point rpt_{tid} are linearly ordered

    // load an allele into a register
    AlleleType allele = alleles[ blockIdx.x * 1024 + tid ];
    __syncthreads();

    // casting to int to allow for negative indices
    int min = 0, max = (eEnd - eStart) - 1;

    // binary search of recombination point list
    while( min <= max ) {
        unsigned int mid = ((max - min) / 2) + min;
        
        if( recomb_points[ mid ] < allele ) {
            min = mid + 1;
        } else if( recomb_points[ mid ] > allele ) {
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

    for( unsigned int i = 1; i < CC::WARP_SIZE; i <<= 1) {
        unsigned int m = __shfl_down(mask, i, CC::WARP_SIZE);
        if( !(threadIdx.x & ((i << 1) - 1)) ) mask |= m;
    }

    if(threadIdx.x == 0 ) {
        sequences[ blockIdx.y * sequence_width + blockIdx.x * blockDim.x + threadIdx.y ] = mask;
    }
}

void crossover::operator()(  real_type         * rand_pool
                            , allele_type       * allele_list
                            , event_count_type  * event_list
                            , int_type          * sequences
                            , size_type nSequences
                            , size_type nAlleles
                            , size_type sequence_width  ) {
    assert( nAlleles % comp_cap_type::THREADS_PER_BLOCK == 0 );
    assert( sequence_width >= nAlleles / comp_cap_type::WARP_SIZE );

    while( nSequences ) {
        unsigned int N = nSequences;
        if( N > MAX_EVENTS ) { N = MAX_EVENTS; }
        cudaMemcpyToSymbol( g_cEvents, event_list, (N + 1) * sizeof( event_count_type ), 0, cudaMemcpyDeviceToDevice );

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

        // assumes the execution configuration specifies sequences to generate
        // allele tile offset determined from execution configuration
        crossover_matrix< real_type, comp_cap_type ><<< blocks, threads >>>( rand_pool, allele_list, sequences, sequence_width );


        std::cerr << "Blocks: < " << blocks.x << ", " << blocks.y << ", " << blocks.z << " >" << std::endl;
        std::cerr << "Threads: < " << threads.x << ", " << threads.y << ", " << threads.z << " >" << std::endl;
        cudaDeviceSynchronize();

        event_list += N;
        nSequences -= N;
        sequences += (N * sequence_width);
    }
}

#endif  // CROSSOVER_MATRIX_5_CUH_
