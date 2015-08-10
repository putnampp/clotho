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
#ifndef CROSSOVER_MATRIX_4_CUH_
#define CROSSOVER_MATRIX_4_CUH_

#include "clotho/cuda/crossover/crossover_matrix_def.hpp"
#include "clotho/cuda/compute_capability.hpp"

#include <cuda.h>

template < >
class crossover< 4 > {
public:
    typedef double                      real_type;
    typedef double                      allele_type;
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


/**
 *  Execution Configuration:
 *      Blocks:
 *          - x => 1 block per THREAD_COUNT alleles in the parent population
 *          - y => 1 block
 *          - z => 1 block
 *      Threads: 1024 threads (compute_capability< 3, 0 >)
 *          - x => 32 threads per warp
 *          - y => 32 warps per block
 *          - z => 1
 **/
template < class RealType, class CC >
__global__ void crossover_matrix_4( RealType * rand_pool
                                    , RealType * allele_list
                                    , unsigned int * evt_list
                                    , unsigned int * sequences
                                    , unsigned int nSequences
                                    , unsigned int nAlleles
                                    , unsigned int sequence_width ) {
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    __shared__ unsigned int sBuffer[ 1024 ];
    __shared__ RealType sRands[ 1024 ];

    sBuffer[tid] = 0;

    RealType allele = allele_list[ blockIdx.x * 1024 + tid ];
    __syncthreads();

    // consider sorting alleles within warps (and thread masks)
    // warp sort is quick (bitonic) and once per sequence block
    // Theory: sorting the alleles results in fewer bank conflicts when
    // accessing event counts from share memory for threads within the
    // same warp.

    unsigned int allele_bin_lo = allele * 1024;
    unsigned int allele_bin_hi = allele_bin_lo + 1;

    sequences += blockIdx.x * blockDim.x;    // block relative shift of sequence array

    while( nSequences-- ) {
        // assumption: this will be a single read and broadcast
        // assumption: evt_list has nSequences + 1 elements

        // assumption: 1 event bin/thread
        // load event count into thread register
        unsigned int eCount = evt_list[ tid ];
        __syncthreads();

        // prefix sum of event count in thread registers
        // sum within warps
#pragma unroll
        for( unsigned int i = 1; i < 32; i <<=1 ) {
            unsigned int t = __shfl_up( eCount, i, 32 );
            eCount += (!!( threadIdx.x >= i )) * t; // double negation should gaurantee that coefficient is the integer value 0 or 1; inline Kronecker delta function
                                                        // attempting to avoid branch divergence; though may not be an issue in this simple case
                                                        // uncertain whether using integer multiplication instead of inline conditional logic is more efficient
        }

        // share partial sum with other warps in block
        if( threadIdx.x == 31 ) {
             sBuffer[ threadIdx.y ] = eCount;
        }
        __syncthreads(); // all warps need to finish sharing 

        // broadcast partial sums to all warps
        unsigned int partial_sum = sBuffer[ threadIdx.x ];

#pragma unroll
        for( unsigned int i = 1; i < 32; i <<= 1 ) {
            unsigned int t = __shfl_up( partial_sum, i, 32 );
            partial_sum += (!!(threadIdx.x >= i )) * t;
        }

        unsigned int t = __shfl( partial_sum, threadIdx.y - 1, 32 );
        eCount += (!!(threadIdx.y != 0)) * t;

        // store prefix sum in shared buffer
        sBuffer[ tid ] = eCount;
        __syncthreads();

        // within each warp, partial_sum of max warp thread is total events
        // use of __shfl to avoid shared memory read
        t = __shfl( partial_sum, 31, 32 );
        partial_sum = t;

        if( partial_sum ) {  // true for all threads in block

            // copy uniform randoms from global memory
            if( tid < partial_sum ) {
                sRands[ tid ] = rand_pool[ tid ];
            }
            __syncthreads();

            rand_pool += partial_sum;
            
            unsigned int min_events = sBuffer[ allele_bin_lo ];
            unsigned int max_events = sBuffer[ allele_bin_hi ];

            unsigned int N = (max_events - min_events);

            RealType curmax = 0.;
            RealType bin_lo = ((RealType) allele_bin_lo) / 1024.;
            while( N ) {
                // walk the events in a bin
                // count all events that preceed the thread's allele
                RealType r = sRands[ min_events ];
                curmax += (log( r ) / ((RealType) N));

                r = bin_lo + (( 1.0 - exp(curmax)) / 1024.);

                min_events += (!!(r < allele));
                --N;
            }
            __syncthreads();
            unsigned int mask = ((min_events & 1) * (1 << threadIdx.x));

            // reduce masks
            for( unsigned int i = 1; i < 32; i <<=1 ) {
                unsigned int m = __shfl_down( mask, i, 32 );
                mask |= ( (!(threadIdx.x & (( i << 1 ) - 1))) * m );
            }

            if( threadIdx.x == 0 ) {
                sequences[ threadIdx.y ] = mask;
            }
        } else if( threadIdx.y == 0 ) {
            sequences[ threadIdx.x ] = 0;
        }

        __syncthreads();

        evt_list += 1024;
        sequences += sequence_width;
    }
}

void crossover< 4 >::operator()(  real_type         * rand_pool
                            , allele_type       * allele_list
                            , event_count_type  * evt_list
                            , int_type          * sequences
                            , size_type nSequences
                            , size_type nAlleles
                            , size_type sequence_width  ) {
    assert( nAlleles % comp_cap_type::THREADS_PER_BLOCK == 0 );
    assert( sequence_width >= nAlleles / comp_cap_type::WARP_SIZE );

    // execution configuration:
    // blockIdx.x = group of THREAD_PER_BLOCK alleles (1 column tile per 1024 alleles)
    // blockIdx.y = offset of offspring sequence (1 - row tile per offspring sequence)
    // blockIdx.z = 1; unused
    //
    // 1 allele/thread therefore
    // threadIdx.x = warp offset
    // threadIdx.y = warps per block
    // threadIdx.z = 1; unused
    dim3 blocks( nAlleles / comp_cap_type::THREADS_PER_BLOCK, 1, 1)
            , threads( comp_cap_type::WARP_SIZE, comp_cap_type::WARP_PER_BLOCK, 1);

    crossover_matrix_4< real_type, comp_cap_type ><<< blocks, threads >>>( rand_pool, allele_list, evt_list, sequences, nSequences, nAlleles, sequence_width );

}
#endif  // CROSSOVER_MATRIX_4_CUH_