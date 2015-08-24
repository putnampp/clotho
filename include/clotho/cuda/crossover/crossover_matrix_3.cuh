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
#ifndef CROSSOVER_MATRIX_3_CUH_
#define CROSSOVER_MATRIX_3_CUH_

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

#include "clotho/cuda/crossover/poisson_distribution.hpp"

template < >
class crossover< 3 > {
public:
    typedef double                      real_type;
    typedef double                      allele_type;
    typedef unsigned int                event_count_type;
    typedef unsigned int               int_type;
    typedef unsigned int                size_type;
    typedef compute_capability< 3, 0 >  comp_cap_type;

    static const unsigned int ALLELE_PER_INT = 32;
    static const unsigned int MAX_EVENTS = ((comp_cap_type::MAX_CONSTANT_MEMORY / sizeof( event_count_type )) >> 1);    // use half of the constant space for event counts

    typedef unsigned long long      seed_type;
    typedef curandStateMtgp32_t     state_type;
    typedef mtgp32_kernel_params_t  state_param_type;

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
    seed_type m_seed;
};

/**
 * The use of texture memory for the alleles does not seem to reduce runtime
 * Preliminary tests indicate that the texture memory adds some over to the process
 */
#ifdef USE_TEXTURE_MEMORY_FOR_ALLELE
texture< float, 1, cudaReadModeElementType > allele_tex;
#endif  // USE_TEXTURE_MEMORY_FOR_ALLELE

//#ifndef USE_PARTIAL_RANDOM_POOL
//#define USE_PARTIAL_RANDOM_POOL
//#endif  // USE_PARTIAL_RANDOM_POOL

template < class RealType >
__global__ void random_test( curandStateMtgp32_t * states, RealType * rands ) {
    int tid = threadIdx.y * blockDim.x + threadIdx.x;

    RealType r = curand_mtgp32_single_specific( &states[blockIdx.x], tid, blockDim.x ) - 1.0;

    rands[ blockIdx.x * blockDim.x + tid ] = r;
}

template < class RealType >
__global__ void random_test( curandStateMtgp32_t * states, RealType * rands, unsigned int N ) {
    int tid = threadIdx.y * blockDim.x + threadIdx.x;

    while( N ) {
        RealType r = curand_mtgp32_single_specific( &states[blockIdx.x], tid, blockDim.x );

        if( tid < N ) {
            rands[ tid ] = r;
        }
        __syncthreads();

        rands += blockDim.x;
        N -= ((N > blockDim.x) ? blockDim.x : N );
    }
}

template < class RealType >
__global__ void crossover_kernel_3( curandStateMtgp32_t * states
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

    unsigned int i;

    __shared__ RealType rand_pool[ THREAD_NUM ];
    __shared__ unsigned int event_hash[ THREAD_NUM + 1];

    rho /= 256.0;

    event_hash[tid] = 0;    // clear event_hash

#ifdef LOG_RANDOM_EVENTS
    pool += blockIdx.x * THREAD_NUM;
    evt_list += blockIdx.x * THREAD_NUM;
#endif  // LOG_RANDOM_EVENTS

#ifdef USE_PARTIAL_RANDOM_POOL
    unsigned int rand_offset = THREAD_NUM, _count;
#endif  // USE_PARTIAL_RANDOM_POOL

    int seq_idx = blockIdx.x;
    while( seq_idx < nSequences ) {

        unsigned int * seq = sequences + (seq_idx * sequence_width);

        event_hash[ tid + 1 ] = 0;  // clear the event_hash
        __syncthreads();

        unsigned int rand = curand_poisson( &states[blockIdx.x], rho );
        __syncthreads();

        for( i = 1; i < 32; i <<= 1 ) {
            unsigned int r = __shfl_up( rand, i );
            rand += ( (lane_id >= i ) * r );
        }

        if( lane_id == 31 ) {
            event_hash[ 32 + warp_id ] = rand;
        }
        __syncthreads();

        unsigned int _sum = event_hash[ 32 + lane_id ];
        for( i = 1; i < (THREAD_NUM >> 5); i <<= 1 ) {
            unsigned int s = __shfl_up( _sum, i );
            _sum += (( lane_id >= i ) * s);
        }

        unsigned int s = __shfl( _sum, warp_id - 1);
        rand += (( warp_id != 0 ) * s);
        __syncthreads();

        event_hash[ tid + 1 ] = rand;

#ifdef LOG_RANDOM_EVENTS
        evt_list[tid] = rand;
        evt_list += BLOCK_NUM_MAX * THREAD_NUM;
#endif  // LOG_RANDOM_EVENTS

        __syncthreads();

        i = event_hash[tid];    // minimum event index
        __syncthreads();

#ifdef USE_PARTIAL_RANDOM_POOL
        _count = __shfl( _sum, 31 ); // total number of events for sequence
        if( rand_offset + _count >= THREAD_NUM ) {
            rand_pool[ tid ] = curand_uniform( &states[blockIdx.x] );
            rand_offset = 0;
        }
        i += rand_offset;
        rand += rand_offset;
#else
        rand_pool[ tid ] = curand_uniform( &states[blockIdx.x] );
#endif  // USE_PARTIAL_RANDOM_POOL
        __syncthreads();

        RealType accum = 0.;
        while (i < rand) {
            RealType t = rand_pool[ i ];

            accum += (log( t ) / (RealType)(rand - i));

            rand_pool[i++] = ((((RealType)tid) + (1.0 - exp(accum))) / 256.0);
        }
        __syncthreads();

#ifdef LOG_RANDOM_EVENTS
#ifdef USE_PARTIAL_RANDOM_POOL
        if( rand_offset <= tid && tid < rand_offset + _count ) {
            pool[tid - rand_offset] = rand_pool[tid];
        }
#else
        pool[tid] = rand_pool[ tid ];
#endif  // USE_PARTIAL_RANDOM_POOL

        pool += BLOCK_NUM_MAX * THREAD_NUM;
        __syncthreads();
#endif  // LOG_RANDOM_EVENTS

        i = tid;
        while( i < nAlleles ) {
#ifdef USE_TEXTURE_MEMORY_FOR_ALLELE
            RealType _allele = tex1Dfetch( allele_tex, i );
#else
            RealType _allele = allele_list[ i ];
#endif  // USE_TEXTURE_MEMORY_FOR_ALLELE

            rand = (unsigned int) (_allele * 256.0);

            unsigned int e_min = event_hash[ rand++ ];
            _sum = event_hash[ rand ];

            unsigned int c = e_min;

#ifdef USE_PARTIAL_RANDOM_POOL
            e_min += rand_offset;
            _sum += rand_offset;
#endif  // USE_PARTIAL_RANDOM_POOL

            while( e_min < _sum ) {
                accum = rand_pool[ e_min++ ];
                c += ( _allele > accum);
            }
            __syncthreads();

            c = ((c & 1) * (1 << lane_id));

            for( rand = 1; rand < 32; rand <<= 1 ) {
                e_min = __shfl_down( c, rand );
                c |= ((!( tid & (( rand << 1) - 1))) * e_min);
            }

            if( lane_id == 0) {
                seq[ warp_id ] = c;
            }
            __syncthreads();

            i += THREAD_NUM;
            seq += (THREAD_NUM >> 5);
        }
        __syncthreads();

        seq_idx += BLOCK_NUM_MAX;
#ifdef USE_PARTIAL_RANDOM_POOL
        rand_offset += _count;
#endif  // USE_PARTIAL_RANDOM_POOL
    }
}

crossover< 3 >::crossover( ) :
    dStates( NULL )
    , dParams( NULL)
{
    initialize();
}

void crossover< 3 >::initialize( ) {
    // defines from curand_mtgp32.h in curand library (as of CUDA 6.5)
    // THREAD_NUM == MTGPDC_FLOOR_2P == 256
    // BLOCK_NUM_MAX == CURAND_NUM_MTGP32_PARAMS == 200

    m_seed = clotho::utility::clock_type::now().time_since_epoch().count();

    assert( cudaMalloc( (void ** ) &dStates, STREAM_COUNT * BLOCK_NUM_MAX * sizeof( state_type ) ) == cudaSuccess );
    assert( cudaMalloc( (void ** ) &dParams, sizeof( state_param_type ) ) == cudaSuccess );

    assert( curandMakeMTGP32Constants( MTGPDC_PARAM_TABLE, dParams ) == CURAND_STATUS_SUCCESS );

    assert( curandMakeMTGP32KernelState( dStates, MTGPDC_PARAM_TABLE, dParams, BLOCK_NUM_MAX, m_seed) == CURAND_STATUS_SUCCESS );
}

void crossover< 3 >::operator()(  real_type * rand_pool
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

    crossover_kernel_3<<< BLOCK_NUM_MAX, THREAD_NUM >>>( dStates, rand_pool, allele_list, event_list, sequences, nSequences, nAlleles, sequence_width, rho );
}

void crossover< 3 >::generate_test( real_type * rand_pool, size_t N ) {
    while( N ) {
        unsigned int T = THREAD_NUM;
        unsigned int B = N / T;
        if( B ) {
            B = (( B < BLOCK_NUM_MAX ) ? B : BLOCK_NUM_MAX );
        } else {
            B = 1;
            T = N;
        }

        random_test<<< B, T >>>( dStates, rand_pool );
        rand_pool += B * T;
        N -= B * T;
    }
}

crossover< 3 >::~crossover() {
    cudaFree( dStates );
    cudaFree( dParams );
}

void crossover< 3 >::get_state( boost::property_tree::ptree & n ) {
    n.put( "crossover.version", 3 );
    n.put( "crossover.curand.state_type", clotho::cuda::curand_helper< typename crossover< 3 >::state_type >::StateName );
    
    n.add_child( "crossover.curand.seed", m_seed  );

#ifdef USE_TEXTURE_MEMORY_FOR_ALLELE
    n.put( "crossover.use_texture_memory_for_allele", true );
#else
    n.put( "crossover.use_texture_memory_for_allele", false );
#endif  // USE_TEXTURE_MEMORY_FOR_ALLELE
}

#endif  // CROSSOVER_MATRIX_3_CUH_
