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

    void adjust_alleles( allele_type * allele_list, size_type N );

    void operator()(  real_type * rand_pool
                    , allele_type       * allele_list
                    , event_count_type  * event_list
                    , int_type          * sequences
                    , size_type nSequences
                    , size_type nAlleles
                    , size_type sequence_width
                    , real_type rho = 0.1 );

    void get_state( boost::property_tree::ptree & s );

    virtual ~crossover();

protected:
    state_type * dStates;
    state_param_type * dParams;
    std::vector< unsigned long long > m_seeds;

    cudaStream_t streams[ STREAM_COUNT ];
};

template < class IntType >
__device__ IntType warp_event_range( IntType count, IntType lane_id, IntType & lo, IntType & hi ) {
    IntType max = count;
    hi = count;

    for( unsigned int i = 1; i < 32; i <<= 1 ) {
        IntType tmp = __shfl_up( max, i );
        max = ((tmp > max) ? tmp : max);

        tmp = __shfl_up( hi , i );
        hi += (( lane_id >= i ) * tmp);
    }

    lo = __shfl_up( hi, 1 );
    lo *= (lane_id != 0);

    max = __shfl( max, 31 );
    return max;
}

template < class IntType >
__global__ void get_warp_event_range( IntType * evt, IntType * out ) {
    IntType tid = threadIdx.y * blockDim.x + threadIdx.x;
    IntType lane_id = (tid & 31);

    IntType max = evt[lane_id];

    IntType evt_lo, evt_hi;
    max = warp_event_range( max, lane_id, evt_lo, evt_hi );
    __syncthreads();

    out[ lane_id ] = evt_lo;
    lane_id += 32;

    out[ lane_id ] = evt_hi;
    lane_id += 32;

    out[ lane_id ] = max;
}

// every warp:
// 1) shuffle up the maximum number of events
// 2) compute the prefix sum of events per bin
template < class RealType >
__global__ void remap_alleles( RealType * allele_list, size_t N ) {
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    RealType bin_start = ((RealType)(tid & 31) / (RealType)32);

    unsigned int a_idx = tid;
    while( a_idx < N ) {
        RealType a = allele_list[ a_idx ];

        a = bin_start + (a / (RealType)32);

        allele_list[ a_idx ] = a;

        a_idx += blockDim.x * blockDim.y;
    }
}

template < class IntType >
__device__ void persist_mask( IntType mask, IntType warp_id, IntType lane_id, IntType * seq ) {
    // collapse masks to single crossover mask per warp
    // mask will exist in lane 0 for all warps
#pragma unroll
    for( unsigned int j = 1; j < 32; j <<= 1 ) {
        unsigned int tmp = __shfl_down( mask, j );
        mask |= ((!( lane_id & (( j << 1) - 1))) * tmp);
    }

    // use single thread per warp to write/store
    // crossover mask to global memory
    if( lane_id == 0) {
        seq[ warp_id ] = mask;
    }
}

template < class IntType >
__device__ void persist_mask( unsigned char mask, IntType warp_id, IntType lane_id, unsigned char * seq ) {
    seq[ warp_id * 32 + lane_id ] = mask;
}

template < class RealType, class IntType >
__device__ IntType count_preceeding_lane_events( volatile RealType * events, RealType allele, IntType lane_id, IntType warp_max, IntType lane_max ) {
    IntType c = 0;
    while( warp_max-- ) {
        RealType e = events[ lane_id ];

        c += (( lane_max ) && (e < allele));

        lane_max -= (!!(lane_max));
        lane_id += 32;
    }

    return c;
}

template < class RealType, class IntType >
__device__ void fill_random_pool( curandStateMtgp32_t * state, volatile RealType * rands, RealType bin_start, RealType bin_width, IntType offset, IntType N ) {
    N += (!!(N & (THREAD_NUM - 1))) * (THREAD_NUM - (N & (THREAD_NUM - 1))); // pad N st N = 256 n; assumes THREAD_NUM = 2^{p}

    // over generate random numbers per bin
    // avoids branch divergence
    while(offset < N ) {
        RealType event = curand_uniform( state );

        rands[ offset ] = bin_start + event * bin_width;

        offset += THREAD_NUM;
    }
}

/**
 *  This version assumes that allele_list is organized into units of WARP_SIZE (32) alleles
 *  that are ordered such that:
 *      allele_unit_idx = floor( allele_location * 32.0 )
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

    typedef unsigned int int_type;

    __shared__ RealType rand_pool[ 1024 ];
    __shared__ int_type event_bins[ THREAD_NUM ];
    int_type tid = threadIdx.y * blockDim.x + threadIdx.x;
    int_type lane_id = (tid & 31);

    int_type i;

    RealType bin_range = (1. / 32.0);
    RealType bin_start =  (RealType) lane_id * bin_range;

    states = &states[blockIdx.x];
    rho *= bin_range; // Warp Size

#ifdef LOG_RANDOM_EVENTS
    pool += blockIdx.x * THREAD_NUM;
    evt_list += blockIdx.x * THREAD_NUM;
#endif  // LOG_RANDOM_EVENTS

    int_type seq_idx = blockIdx.x;
    while( seq_idx < nSequences ) {

        unsigned int * seq = sequences + (seq_idx * sequence_width);

        // use all threads to generate random numbers
        event_bins[ tid ] = curand_poisson( states, rho );
        __syncthreads();

#ifdef LOG_RANDOM_EVENTS
        evt_list[ tid ] = event_bins[ tid ];
        evt_list += BLOCK_NUM_MAX * THREAD_NUM;
#endif  // LOG_RANDOM_EVENTS

        int_type evt_lo, evt_hi;
        int_type max_rounds = warp_event_range( event_bins[lane_id], lane_id, evt_lo, evt_hi);
        __syncthreads();

        // within a warp each thread (bin) has [evt_lo, evt_hi) events

        // over generate random numbers per bin
        // avoids branch divergence
        fill_random_pool( states, rand_pool, bin_start, bin_range, lane_id, (max_rounds << 5) );
        __syncthreads();

#ifdef LOG_RANDOM_EVENTS
        i = tid;
        while( i < 1024 ) {
            pool[ i ] = rand_pool[ i ];
            i += THREAD_NUM;
        }
        pool += BLOCK_NUM_MAX * 1024;
        __syncthreads();
#endif  // LOG_RANDOM_EVENTS

        i = tid;
        while( i < nAlleles ) {
            RealType _allele = allele_list[ i ];

            unsigned int cmask = count_preceeding_lane_events( rand_pool, _allele, lane_id, max_rounds, (evt_hi - evt_lo) );
            __syncthreads();

            cmask = (((cmask + evt_lo) & 1) * (1 << lane_id));  // translate event count to lane relative crossover mask

///////////////////////////////////////////////////////////
            persist_mask( cmask, (tid >> 5), lane_id, seq );
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

void crossover< 2 >::adjust_alleles( allele_type * allele_list, size_type N ) {
    remap_alleles<<< 1, 1024 >>>( allele_list, N );
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

void crossover< 2 >::get_state( boost::property_tree::ptree & n ) {
    n.put( "crossover.version", 2 );
    n.put( "crossover.curand.state_type", clotho::cuda::curand_helper< typename crossover< 2 >::state_type >::StateName );
    
    boost::property_tree::ptree sds;
    clotho::utility::add_value_array( sds, m_seeds.begin(), m_seeds.end() );

    n.add_child( "crossover.curand.seed", sds  );

#ifdef USE_TEXTURE_MEMORY_FOR_ALLELE
    n.put( "crossover.use_texture_memory_for_allele", true );
#else
    n.put( "crossover.use_texture_memory_for_allele", false );
#endif  // USE_TEXTURE_MEMORY_FOR_ALLELE
}

crossover< 2 >::~crossover() {
    cudaFree( dStates );
    cudaFree( dParams );

    for( unsigned int i = 0; i < STREAM_COUNT; ++i ) {
        cudaStreamSynchronize( streams[i] );

        cudaStreamDestroy( streams[i] );
    }
}

#endif  // CROSSOVER_MATRIX_2_CUH_
