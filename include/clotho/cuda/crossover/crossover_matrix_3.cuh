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
//#include <curand_poisson.h>
#include <curand_uniform.h>

#include <cuda_runtime.h>

#include "clotho/cuda/curand_helper.hpp"

#include "clotho/utility/timer.hpp"
#include "clotho/utility/log_helper.hpp"

#include "clotho/cuda/crossover/poisson_distribution.hpp"
#include "clotho/cuda/crossover/persist_sequence.hpp"

template < >
class crossover< 3 > {
public:
    typedef float                      real_type;
    typedef float                      allele_type;
    typedef unsigned int                event_count_type;
    typedef unsigned int               int_type;
    typedef unsigned int                size_type;
    typedef compute_capability< 3, 0 >  comp_cap_type;

    static const unsigned int ALLELE_PER_INT = 32;
    static const unsigned int MAX_EVENTS = ((comp_cap_type::MAX_CONSTANT_MEMORY / sizeof( event_count_type )) >> 1);    // use half of the constant space for event counts

    typedef unsigned long long      seed_type;

#ifdef USE_MERSENNE_TWISTER
    typedef curandStateMtgp32_t     state_type;
    typedef mtgp32_kernel_params_t  state_param_type;

    static const int_type   BLOCKS_PER_KERNEL = BLOCK_NUM_MAX;
    static const int_type   THREADS_PER_BLOCK = THREAD_NUM;
#else
    typedef curandState_t           state_type;

    static const int_type   BLOCKS_PER_KERNEL = 200;
    static const int_type   THREADS_PER_BLOCK = 128;
#endif  // USE_MERSENNE_TWISTER

    typedef poisson_cdf< real_type, 32 > poisson_type;

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

#ifdef USE_MERSENNE_TWISTER
    state_param_type * dParams;
#endif  // USE_MERSENNE_TWISTER

    seed_type m_seed;

    poisson_type * m_pois_dist;
};

/**
 * The use of texture memory for the alleles does not seem to reduce runtime
 * Preliminary tests indicate that the texture memory adds some over to the process
 */
#ifdef USE_TEXTURE_MEMORY_FOR_ALLELE
texture< float, 1, cudaReadModeElementType > allele_tex;
#endif  // USE_TEXTURE_MEMORY_FOR_ALLELE

template < class StateType, class RealType >
__global__ void random_test( StateType * states, RealType * rands ) {
    int tid = threadIdx.y * blockDim.x + threadIdx.x;

    RealType r = curand_uniform( &states[blockIdx.x * blockDim.x + tid] );

    rands[ blockIdx.x * blockDim.x + tid ] = r;
}

template < class StateType, class RealType >
__global__ void random_test( StateType * states, RealType * rands, unsigned int N ) {
    int tid = threadIdx.y * blockDim.x + threadIdx.x;

    while( N ) {
        RealType r = curand_uniform( &states[blockIdx.x * blockDim.x + tid] );

        if( tid < N ) {
            rands[ tid ] = r;
        }
        __syncthreads();

        rands += blockDim.x;
        N -= ((N > blockDim.x) ? blockDim.x : N );
    }
}

template < class real_type >
__global__ void random_test( curandStateMtgp32_t * states, real_type * rands ) {
    int tid = threadIdx.y * blockDim.x + threadIdx.x;

    real_type r = curand_mtgp32_single_specific( &states[blockIdx.x], tid, blockDim.x ) - 1.0;

    rands[ blockIdx.x * blockDim.x + tid ] = r;
}

template < class real_type >
__global__ void random_test( curandStateMtgp32_t * states, real_type * rands, unsigned int N ) {
    int tid = threadIdx.y * blockDim.x + threadIdx.x;

    while( N ) {
        real_type r = curand_mtgp32_single_specific( &states[blockIdx.x], tid, blockDim.x );

        if( tid < N ) {
            rands[ tid ] = r;
        }
        __syncthreads();

        rands += blockDim.x;
        N -= ((N > blockDim.x) ? blockDim.x : N );
    }
}

template < class XOVER >
__global__ void crossover_kernel_3( typename XOVER::state_type * states
                                    , typename XOVER::real_type * pool
                                    , typename XOVER::real_type * allele_list
                                    , typename XOVER::int_type * evt_list
                                    , typename XOVER::int_type * sequences
                                    , typename XOVER::int_type nSequences
                                    , typename XOVER::int_type nAlleles
                                    , typename XOVER::int_type sequence_width
                                    , typename XOVER::poisson_type * pois ) {
    if( gridDim.x != XOVER::BLOCKS_PER_KERNEL || blockDim.x != XOVER::THREADS_PER_BLOCK ) return;

    typedef typename XOVER::real_type   real_type;
    typedef typename XOVER::state_type  state_type;
    typedef typename XOVER::int_type    int_type;

    int_type tid = threadIdx.y * blockDim.x + threadIdx.x;
    int_type lane_id = (tid & 31);
    int_type warp_id = (tid >> 5);

    int_type i;

    __shared__ real_type    s_pois_cdf[ XOVER::poisson_type::MAX_K ];
    __shared__ real_type    rand_pool[ XOVER::THREADS_PER_BLOCK ];
    __shared__ unsigned int event_hash[ XOVER::THREADS_PER_BLOCK + 1];

    event_hash[ tid ] = 0;

    if( tid < XOVER::poisson_type::MAX_K ) {
        s_pois_cdf[ tid ] = pois->_cdf[tid];
    }
    unsigned int max_k = pois->max_k;
    __syncthreads();

#ifndef USE_MERSENNE_TWISTER
    typename XOVER::state_type local_state = states[ blockIdx.x * blockDim.x *blockDim.y + tid ];
#endif  // USE_MERSENNE_TWISTER

#ifdef LOG_RANDOM_EVENTS
    pool += blockIdx.x * XOVER::THREADS_PER_BLOCK;
    evt_list += blockIdx.x * XOVER::THREADS_PER_BLOCK;
#endif  // LOG_RANDOM_EVENTS

    int seq_idx = blockIdx.x;
    while( seq_idx < nSequences ) {

        unsigned int * seq = sequences + (seq_idx * sequence_width);

#ifdef USE_MERSENNE_TWISTER
        real_type x = curand_uniform( &states[blockIdx.x] );
        rand_pool[ tid ] = curand_uniform( &states[blockIdx.x] );
#else
        real_type x = curand_uniform( &local_state );
        rand_pool[ tid ] = curand_uniform( &local_state );
#endif  // USE_MERSENNE_TWISTER

        unsigned int rand = _find_poisson_maxk32( s_pois_cdf, x, max_k );
        __syncthreads();

#ifdef LOG_RANDOM_EVENTS
        evt_list[tid] = rand;
        evt_list += XOVER::BLOCKS_PER_KERNEL * XOVER::THREADS_PER_BLOCK;
#endif  // LOG_RANDOM_EVENTS

        for( i = 1; i < 32; i <<= 1 ) {
            unsigned int r = __shfl_up( rand, i );
            rand += ( (lane_id >= i ) * r );
        }

        if( lane_id == 31 ) {
            event_hash[ 32 + warp_id ] = rand;
        }
        __syncthreads();

        unsigned int _sum = event_hash[ 32 + lane_id ];
        for( i = 1; i < (XOVER::THREADS_PER_BLOCK >> 5); i <<= 1 ) {
            unsigned int s = __shfl_up( _sum, i );
            _sum += (( lane_id >= i ) * s);
        }

        unsigned int s = __shfl( _sum, warp_id - 1);
        rand += (( warp_id != 0 ) * s);
        __syncthreads();

        event_hash[ tid + 1 ] = rand;
        __syncthreads();

        i = event_hash[tid];    // minimum event index
        __syncthreads();

        real_type accum = 0.;
        while (i < rand) {
            real_type t = rand_pool[ i ];

            accum += (log( t ) / (real_type)(rand - i));

            rand_pool[i++] = ((((real_type)tid) + (1.0 - exp(accum))) / ((real_type)XOVER::THREADS_PER_BLOCK));
        }
        __syncthreads();

#ifdef LOG_RANDOM_EVENTS
        pool[tid] = rand_pool[ tid ];
        pool += XOVER::BLOCKS_PER_KERNEL * XOVER::THREADS_PER_BLOCK;
        __syncthreads();
#endif  // LOG_RANDOM_EVENTS

        i = tid;
        while( i < nAlleles ) {
#ifdef USE_TEXTURE_MEMORY_FOR_ALLELE
            real_type _allele = tex1Dfetch( allele_tex, i );
#else
            real_type _allele = allele_list[ i ];
#endif  // USE_TEXTURE_MEMORY_FOR_ALLELE

            rand = (unsigned int) (_allele * ((real_type)XOVER::THREADS_PER_BLOCK));

            unsigned int e_min = event_hash[ rand++ ];
            _sum = event_hash[ rand ];

            int_type cmask = e_min;

            while( e_min < _sum ) {
                accum = rand_pool[ e_min++ ];
                cmask += ( _allele > accum);
            }
            __syncthreads();

            cmask = ((cmask & 1) * (1 << lane_id));

            persist_mask_unrolled( cmask, warp_id, lane_id, seq );
            __syncthreads();

            i += XOVER::THREADS_PER_BLOCK;
            seq += (XOVER::THREADS_PER_BLOCK >> 5);
        }
        __syncthreads();

        seq_idx += XOVER::BLOCKS_PER_KERNEL;
    }
#ifndef USE_MERSENNE_TWISTER
    states[ blockIdx.x * blockDim.x *blockDim.y + tid ] = local_state;
#endif  // USE_MERSENNE_TWISTER
}

crossover< 3 >::crossover( ) :
    dStates( NULL )
#ifdef USE_MERSENNE_TWISTER
    , dParams( NULL)
#endif  // USE_MERSENNE_TWISTER
{
    initialize();
}

void crossover< 3 >::initialize( ) {
    // defines from curand_mtgp32.h in curand library (as of CUDA 6.5)
    // THREAD_NUM == MTGPDC_FLOOR_2P == 256
    // BLOCK_NUM_MAX == CURAND_NUM_MTGP32_PARAMS == 200

    m_seed = clotho::utility::clock_type::now().time_since_epoch().count();

    assert( cudaMalloc( (void ** ) &m_pois_dist, sizeof( poisson_type ) ) == cudaSuccess );

#ifdef USE_MERSENNE_TWISTER
    assert( cudaMalloc( (void ** ) &dStates, BLOCKS_PER_KERNEL * sizeof( state_type ) ) == cudaSuccess );
    assert( cudaMalloc( (void ** ) &dParams, sizeof( state_param_type ) ) == cudaSuccess );

    assert( curandMakeMTGP32Constants( MTGPDC_PARAM_TABLE, dParams ) == CURAND_STATUS_SUCCESS );

    assert( curandMakeMTGP32KernelState( dStates, MTGPDC_PARAM_TABLE, dParams, BLOCKS_PER_KERNEL, m_seed) == CURAND_STATUS_SUCCESS );
#else

    assert( cudaMalloc( (void ** ) &dStates, BLOCKS_PER_KERNEL * THREADS_PER_BLOCK * sizeof( state_type )) == cudaSuccess );

    clotho::cuda::setup_state_kernel<<< BLOCKS_PER_KERNEL, THREADS_PER_BLOCK >>>( dStates, m_seed );
#endif  // USE_MERSENNE_TWISTER
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

    real_type lambda = rho / ((real_type)THREADS_PER_BLOCK);
    make_poisson_cdf_maxk32<<< 1, 32 >>>( m_pois_dist, lambda );
    crossover_kernel_3< crossover<3> ><<< BLOCKS_PER_KERNEL, THREADS_PER_BLOCK >>>( dStates, rand_pool, allele_list, event_list, sequences, nSequences, nAlleles, sequence_width, m_pois_dist );
}

void crossover< 3 >::generate_test( real_type * rand_pool, size_t N ) {
    while( N ) {
        unsigned int T = THREADS_PER_BLOCK;
        unsigned int B = N / T;
        if( B ) {
            B = (( B < BLOCKS_PER_KERNEL ) ? B : BLOCKS_PER_KERNEL );
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

#ifdef USE_MERSENNE_TWISTER
    cudaFree( dParams );
#endif  // USE_MERSENNE_TWISTER

    cudaFree( m_pois_dist );
}

void crossover< 3 >::get_state( boost::property_tree::ptree & n ) {
    n.put( "crossover.version", 3 );
    n.put( "crossover.curand.state_type", clotho::cuda::curand_helper< typename crossover< 3 >::state_type >::StateName );
    
    n.put( "crossover.curand.seed", m_seed  );

#ifdef USE_TEXTURE_MEMORY_FOR_ALLELE
    n.put( "crossover.use_texture_memory_for_allele", true );
#else
    n.put( "crossover.use_texture_memory_for_allele", false );
#endif  // USE_TEXTURE_MEMORY_FOR_ALLELE
}

#endif  // CROSSOVER_MATRIX_3_CUH_
