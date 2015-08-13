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

template < >
class crossover< 3 > {
public:
    typedef double                      real_type;
    typedef double                      allele_type;
    typedef unsigned int                event_count_type;
    typedef unsigned int               int_type;
    typedef unsigned int                size_type;
    typedef compute_capability< 3, 0 >  comp_cap_type;

    static const unsigned int ALLELE_PER_INT = 1;
    static const unsigned int MAX_EVENTS = ((comp_cap_type::MAX_CONSTANT_MEMORY / sizeof( event_count_type )) >> 1);    // use half of the constant space for event counts

    static const unsigned int STREAM_COUNT = 8;

    typedef curandStateMtgp32_t state_type;
    typedef mtgp32_kernel_params_t state_param_type;

    crossover();

    void initialize( );

    void operator()(  real_type * rand_pool
                    , allele_type       * allele_list
                    , event_count_type  * event_list
                    , int_type          * sequences
                    , size_type nSequences
                    , size_type nAlleles
                    , size_type sequence_width );

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
 * The use of texture memory for the alleles does not seem to reduce runtime
 * Preliminary tests indicate that the texture memory adds some over to the process
 */
#ifdef USE_TEXTURE_MEMORY_FOR_ALLELE
texture< float, 1, cudaReadModeElementType > allele_tex;
#endif  // USE_TEXTURE_MEMORY_FOR_ALLELE

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
                                    , RealType * allele_list
                                    , unsigned int * evt_list
                                    , unsigned int * sequences
                                    , unsigned int nSequences
                                    , unsigned int nAlleles
                                    , unsigned int sequence_width ) {
    if( gridDim.x != BLOCK_NUM_MAX || blockDim.x != THREAD_NUM ) return;

    int tid = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int i;

    __shared__ RealType rand_pool[ THREAD_NUM ];
    __shared__ unsigned int event_hash[ THREAD_NUM + 1];

    int seq_idx = blockIdx.x;
    while( seq_idx < nSequences ) {

        unsigned int * seq = sequences + (seq_idx * sequence_width);

        event_hash[ tid + 1 ] = 0;  // clear the event_hash
        __syncthreads();

        rand_pool[ tid ] = curand_uniform( &states[blockIdx.x] );
        __syncthreads();

        unsigned int rand = curand_poisson( &states[blockIdx.x], 30.0 / 256.0 );

        for( i = 1; i < 32; i <<= 1 ) {
            unsigned int r = __shfl_up( rand, i );
            rand += ( ((tid & 31) >= i ) * r );
        }

        if( (tid & 31) == 31 ) {
            event_hash[ 32 + (tid >> 5) ] = rand;
        }
        __syncthreads();

        unsigned int _sum = event_hash[ 32 + (tid & 31) ];
        for( i = 1; i < 32; i <<= 1 ) {
            unsigned int s = __shfl_up( _sum, i );
            _sum += (( (tid & 31) >= i ) * s);
        }

        i = __shfl( _sum, (tid >> 5) - 1);
        rand += (( tid >= 32 ) * i);

        event_hash[ tid + 1 ] = rand;
        __syncthreads();

        int all_idx = tid;
        while( all_idx < nAlleles ) {
#ifdef USE_TEXTURE_MEMORY_FOR_ALLELE
            RealType _allele = tex1Dfetch( allele_tex, all_idx );
#else
            RealType _allele = allele_list[all_idx];
#endif  // USE_TEXTURE_MEMORY_FOR_ALLELE

            i = (unsigned int) (_allele * 256.0);

            unsigned int e_min = event_hash[ i++ ];
            unsigned int e_max = event_hash[ i ];

            RealType accum = 0.;
            unsigned int c = e_min;
            while( e_min < e_max ) {
                RealType t = rand_pool[ e_min ];

                accum += (log( t ) / (RealType)(e_max - e_min));

                t = ((((RealType)i) + (1.0 - exp(accum))) / 256.0);

                c += (( _allele > t ) * 1);
                ++e_min;
            }
            __syncthreads();    // sync to coalesce global memory write

            seq[tid] = (c & 1);

            all_idx += THREAD_NUM;
            seq += THREAD_NUM;
        }

        seq_idx += BLOCK_NUM_MAX;
    }
}

crossover< 3 >::crossover() :
    dStates( NULL )
    , dParams( NULL)
{
    initialize();
}

void crossover< 3 >::initialize( ) {
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

void crossover< 3 >::operator()(  real_type * rand_pool
                    , allele_type       * allele_list
                    , event_count_type  * event_list
                    , int_type          * sequences
                    , size_type nSequences
                    , size_type nAlleles
                    , size_type sequence_width ) {

#ifdef USE_TEXTURE_MEMORY_FOR_ALLELE
    cudaBindTexture(0, allele_tex, allele_list, nAlleles * 4);
#endif  // USE_TEXTURE_MEMORY_FOR_ALLELE

    unsigned int seqs = (nSequences / STREAM_COUNT) + 1;
    unsigned int i = 0;
    while( nSequences ) {
        unsigned int N = ((nSequences > seqs) ? seqs : nSequences);
        
        crossover_kernel_3<<< BLOCK_NUM_MAX, THREAD_NUM, 0, streams[ i ] >>>( dStates + i * BLOCK_NUM_MAX, allele_list, event_list, dSequences, N, nAlleles, sequence_width );

        sequences += N * sequence_width;
        nSequences -= N;
        ++i;
    }

//    crossover_kernel_3<<< BLOCK_NUM_MAX, THREAD_NUM >>>( dStates, allele_list, event_list, sequences, nSequences, nAlleles, sequence_width );
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

    for( unsigned int i = 0; i < STREAM_COUNT; ++i ) {
        cudaStreamSynchronize( streams[i] );

        cudaStreamDestroy( streams[i] );
    }
}

void crossover< 3 >::get_state( boost::property_tree::ptree & n ) {
    n.put( "crossover.version", 3 );
    n.put( "crossover.curand.state_type", clotho::cuda::curand_helper< typename crossover< 3 >::state_type >::StateName );
    
    boost::property_tree::ptree sds;
    clotho::utility::add_value_array( sds, m_seeds.begin(), m_seeds.end() );

    n.add_child( "crossover.curand.seed", sds  );

#ifdef USE_TEXTURE_MEMORY_FOR_ALLELE
    n.put( "crossover.use_texture_memory_for_allele", true );
#else
    n.put( "crossover.use_texture_memory_for_allele", false );
#endif  // USE_TEXTURE_MEMORY_FOR_ALLELE
}

#endif  // CROSSOVER_MATRIX_3_CUH_
