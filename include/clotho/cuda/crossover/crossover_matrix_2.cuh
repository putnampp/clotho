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
#ifndef CROSSOVER_MATRIX_2_CUH_
#define CROSSOVER_MATRIX_2_CUH_

#include "clotho/cuda/crossover/crossover_matrix_def.hpp"

#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>  // needed to define curandStateXORWOW_t for curand_poisson
#include <curand_mtgp32.h>
#include <curand_mtgp32_host.h>
#include <curand_mtgp32_kernel.h>
#include <curand_mtgp32dc_p_11213.h>
#include <curand_poisson.h>
#include <curand_uniform.h>

#include <cuda_runtime.h>

#include "clotho/cuda/curand_helper.hpp"

#include "clotho/utility/timer.hpp"
#include "clotho/utility/log_helper.hpp"

template < >
class crossover< 2 > {
public:
    typedef double                      real_type;
    typedef double                      allele_type;
    typedef unsigned int                event_count_type;
    typedef unsigned int               int_type;
    typedef unsigned int                size_type;
    typedef compute_capability< 3, 0 >  comp_cap_type;

    static const unsigned int ALLELE_PER_INT = 32;
    static const unsigned int MAX_EVENTS = ((comp_cap_type::MAX_CONSTANT_MEMORY / sizeof( event_count_type )) >> 1);    // use half of the constant space for event counts

    static const unsigned int STREAM_COUNT = 8;

    typedef curandStateMtgp32_t state_type;
    typedef mtgp32_kernel_params_t state_param_type;

    crossover( );

    void initialize( );

    void operator()(  real_type * rand_pool
                    , allele_type       * allele_list
                    , event_count_type  * event_list
                    , int_type          * sequences
                    , size_type nSequences
                    , size_type nAlleles
                    , size_type sequence_width
                    , real_type rho = 0.1 );

    void generate_test( real_type * rand_pool, size_t N );

    void get_state( boost::property_tree::ptree & s );

    virtual ~crossover();

protected:
    state_type * dStates;
    state_param_type * dParams;
    std::vector< unsigned long long > m_seeds;

    cudaStream_t streams[ STREAM_COUNT ];
};


/**
 *  This version assumes that allele_list is organized into units of WARP_SIZE (32) alleles
 *  that are ordered such that:
 *      allele_unit_idx = floor( allele_location / 32.0 )
 */
template < class RealType >
__global__ void crossover_kernel_2( curandStateMtgp32_t * states
                                    , RealType * pool
                                    , RealType * allele_list
                                    , unsigned int * evt_list
                                    , unsigned int * sequences
                                    , unsigned int nSequences
                                    , unsigned int nAlleles
                                    , unsigned int sequence_width
                                    , RealType rho ) {
    if( gridDim.x != BLOCK_NUM_MAX || blockDim.x != THREAD_NUM ) return;

    int tid = threadIdx.y * blockDim.x + threadIdx.x;
    int lane_id = (tid & 31);
    int warp_id = (tid >> 5);
    unsigned int lane_mask = (1 << lane_id);

    unsigned int i;

    RealType bin_start =  ((RealType) lane_id / (RealType) 32);

    __shared__ RealType rand_pool[ 1024 ];
    __shared__ unsigned int event_bins[ THREAD_NUM ];

    rho /= ((RealType) 32); // Warp Size

    event_bins[tid] = 0;    // clear event_bins

    int seq_idx = blockIdx.x;
    while( seq_idx < nSequences ) {

        unsigned int * seq = sequences + (seq_idx * sequence_width);

        // use all threads to generate random numbers
        event_bins[ tid ] = curand_poisson( &states[blockIdx.x], rho );
        __syncthreads();

        unsigned int max_rounds = event_bins[ lane_id ];
        unsigned int psum = max_rounds;
        __syncthreads();

        // every warp:
        // 1) shuffle up the maximum number of events
        // 2) compute the prefix sum of events per bin
        for( i = 1; i < 32; i <<= 1 ) {
            unsigned int tmp = __shfl_up( max_rounds, i );
            max_rounds = ((tmp > max_rounds) ? m : max_rounds);

            tmp = __shfl_up( psum, i );
            psum += (( lane_id >= i ) * tmp);
        }

        max_rounds = __shfl( max_rounds, 31 );  // broadcast the max_event

        unsigned int min_events = __shfl_up( psum, 1 );
        min_events *= (lane_id != 0 );
        __syncthreads();

        // within a warp each thread (bin) has [min_events, psum) events

        i = tid;

        // over generate random numbers per bin
        // avoids branch divergence
        unsigned int rcount = (((max_rounds * 32) / THREAD_NUM) + 1) * THREAD_NUM;
        do {
            RealType event = curand_uniform( &states[blockIdx.x] );

            rand_pool[ i ] = bin_start + (event / (RealType) 32);
            i += THREAD_NUM;
        } while( i < rcount );
        __syncthreads();

        i = tid;
        while( i < nAlleles ) {
            RealType _allele = allele_list[ i ];

            unsigned int cmask = min_events;    // any allele in this (thread) lane will have AT LEAST 'min_events' crossover events preceeding it
            unsigned int j = 0, offset = lane_id;

            // not all bins need max_rounds
            // but early termination would introduce branch divergence
            while( j < max_rounds ) {
                RealType e = rand_pool[ offset ];   // load the next event from the bin

                cmask += ((min_events + j < psum ) && ( e < _allele ));

                ++j;
                offset += 32;
            }
            __syncthreads();

            cmask = ((cmask & 1) * lane_mask);  // translate event count to lane relative crossover mask

///////////////////////////////////////////////////////////
            // collapse masks to single crossover mask per warp
            // mask will exist in lane 0 for all warps
            for( rand = 1; rand < 32; rand <<= 1 ) {
                unsigned int tmp = __shfl_down( cmask, rand );
                cmask |= ((!( tid & (( rand << 1) - 1))) * tmp);
            }

            // use single thread per warp to write/store
            // crossover mask to global memory
            if( lane_id == 0) {
                seq[ warp_id ] = cmask;
            }
            __syncthreads();
///////////////////////////////////////////////////////////

            i += THREAD_NUM;
            seq += (THREAD_NUM >> 5);
        }
        __syncthreads();

        seq_idx += BLOCK_NUM_MAX;
    }
}

crossover< 2 >::crossover( ) :
    dStates( NULL )
    , dParams( NULL)
{
    initialize();
}

void crossover< 2 >::initialize( ) {
    // defines from curand_mtgp32.h in curand library (as of CUDA 6.5)
    // THREAD_NUM == MTGPDC_FLOOR_2P == 256
    // BLOCK_NUM_MAX == CURAND_NUM_MTGP32_PARAMS == 200

    //m_seed = clotho::utility::clock_type::now().time_since_epoch().count();

    m_seeds.reserve( STREAM_COUNT );

    assert( cudaMalloc( (void ** ) &dStates, STREAM_COUNT * BLOCK_NUM_MAX * sizeof( state_type ) ) == cudaSuccess );
    assert( cudaMalloc( (void ** ) &dParams, sizeof( state_param_type ) ) == cudaSuccess );

    assert( curandMakeMTGP32Constants( MTGPDC_PARAM_TABLE, dParams ) == CURAND_STATUS_SUCCESS );


    for( unsigned int i = 0; i < STREAM_COUNT; ++i ) {
        cudaStreamCreate( &streams[i] );

        unsigned long long seed = clotho::utility::clock_type::now().time_since_epoch().count();

        assert( curandMakeMTGP32KernelState( dStates + i * BLOCK_NUM_MAX, MTGPDC_PARAM_TABLE, dParams, BLOCK_NUM_MAX, seed) == CURAND_STATUS_SUCCESS );

        m_seeds.push_back(seed);
    }
}

void crossover< 2 >::operator()(  real_type * rand_pool
                    , allele_type       * allele_list
                    , event_count_type  * event_list
                    , int_type          * sequences
                    , size_type nSequences
                    , size_type nAlleles
                    , size_type sequence_width
                    , real_type rho ) {

#ifdef USE_TEXTURE_MEMORY_FOR_ALLELE
    cudaBindTexture(0, allele_tex, allele_list, nAlleles * 4);
#endif  // USE_TEXTURE_MEMORY_FOR_ALLELE

//    unsigned int seqs = (nSequences / STREAM_COUNT) + 1;
//    unsigned int i = 0;
//    while( nSequences ) {
//        unsigned int N = ((nSequences > seqs) ? seqs : nSequences);
//        
//        crossover_kernel_2<<< BLOCK_NUM_MAX, THREAD_NUM, 0, streams[ i ] >>>( dStates + i * BLOCK_NUM_MAX, allele_list, event_list, dSequences, N, nAlleles, sequence_width );
//
//        sequences += N * sequence_width;
//        nSequences -= N;
//        ++i;
//    }

    crossover_kernel_2<<< BLOCK_NUM_MAX, THREAD_NUM >>>( dStates, rand_pool, allele_list, event_list, sequences, nSequences, nAlleles, sequence_width, rho );
}

#endif  // CROSSOVER_MATRIX_2_CUH_
