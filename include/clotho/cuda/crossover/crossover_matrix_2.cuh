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
#include <thrust/device_vector.h>

#include "clotho/cuda/curand_helper.hpp"

#include "clotho/utility/timer.hpp"
#include "clotho/utility/log_helper.hpp"

#include "clotho/cuda/crossover/poisson_distribution.hpp"

//#ifndef USE_MERSENNE_TWISTER
//#define USE_MERSENNE_TWISTER
//#endif  // USE_MERSENNE_TWISTER

template < >
class crossover< 2 > {
public:
    typedef float                      real_type;
    typedef float                      allele_type;
    typedef unsigned int                event_count_type;
    typedef unsigned int               int_type;
    typedef unsigned int                size_type;
    typedef compute_capability< 3, 0 >  comp_cap_type;

    typedef poisson_cdf< real_type, 32> poisson_type;

    static const unsigned int ALLELE_PER_INT = 32;
    static const unsigned int MAX_EVENTS = ((comp_cap_type::MAX_CONSTANT_MEMORY / sizeof( event_count_type )) >> 1);    // use half of the constant space for event counts

#ifdef USE_MERSENNE_TWISTER
    typedef curandStateMtgp32_t     state_type;
    typedef mtgp32_kernel_params_t  state_param_type;

    static const unsigned int BLOCKS_PER_KERNEL = BLOCK_NUM_MAX;
    static const unsigned int THREADS_PER_BLOCK = THREAD_NUM;
#else
    typedef curandState_t           state_type;

    static const unsigned int BLOCKS_PER_KERNEL = 200;
    static const unsigned int THREADS_PER_BLOCK = 128;
#endif  // USE_MERSENNE_TWISTER

    typedef unsigned long long      seed_type;

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

#ifdef USE_MERSENNE_TWISTER
    state_param_type * dParams;
#endif  // USE_MERSENNE_TWISTER

    seed_type m_seed;

    poisson_type * dPoisCDF;
};

template < class StateType >
__global__ void setup_state_kernel( StateType * states, unsigned long long seed ) {
    int id = blockIdx.x * blockDim.x * blockDim.y + threadIdx.y * blockDim.x + threadIdx.x;

    curand_init( seed, id, 0, &states[id] );
}

// every warp:
// 1) shuffle up the maximum number of events
// 2) compute the prefix sum of events per bin
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
__device__ void persist_mask_old( IntType mask, IntType warp_id, IntType lane_id, IntType * seq ) {
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
__device__ void persist_mask( IntType mask, IntType warp_id, IntType lane_id, IntType * seq ) {
    // collapse masks to single crossover mask per warp
    // mask will exist in lane 0 for all warps

    unsigned int tmp = __shfl_down( mask, 1 );
    mask |= ((!( lane_id & 1)) * tmp);
//    mask |= ((( lane_id & 1) == 0) ? tmp : 0 );

    tmp = __shfl_down( mask, 2);
    mask |= ((!( lane_id & 3)) * tmp);
//    mask |= ((( lane_id & 3) == 0) ? tmp : 0 );
    
    tmp = __shfl_down( mask, 4);
    mask |= ((!( lane_id & 7)) * tmp);
//    mask |= ((( lane_id & 7) == 0) ? tmp : 0 );

    tmp = __shfl_down( mask, 8);
    mask |= ((!( lane_id & 15)) * tmp);
//    mask |= ((( lane_id & 15) == 0) ? tmp : 0 );

    tmp = __shfl_down( mask, 16);
    mask |= ((!( lane_id & 31)) * tmp);
//    mask |= ((( lane_id & 31) == 0) ? tmp : 0 );

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
__device__ IntType count_preceeding_lane_events( volatile RealType * events, RealType allele, IntType lane_id, IntType warp_max ) {
    IntType c = 0;
    while( warp_max-- ) {
        RealType e = events[ lane_id ];

        //c += (e < allele);
        c += ((e < allele) ? 1 : 0);

        lane_id += 32;
    }

    return c;
}

template < class XOVER >
__device__ void fill_random_pool( typename XOVER::state_type * state, 
                                  volatile typename XOVER::real_type * rands, 
                                  typename XOVER::real_type bin_start,
                                  typename XOVER::real_type bin_width,
                                  typename XOVER::int_type offset,
                                  typename XOVER::int_type N,
                                  typename XOVER::int_type M ) {
//    N += (!!(N & (THREAD_NUM - 1))) * (THREAD_NUM - (N & (THREAD_NUM - 1))); // pad N st N = 256 n; assumes THREAD_NUM = 2^{p}

    typename XOVER::int_type tail = (N & (XOVER::THREADS_PER_BLOCK - 1));
    N += ( (!!(tail)) * (XOVER::THREADS_PER_BLOCK - tail));

    // over generate random numbers per bin
    // avoids branch divergence
    while(offset < N ) {
        typename XOVER::real_type event = curand_uniform( state );

        rands[ offset ] = ((M > 0) ? (bin_start + event * bin_width) : 1.0);

        M -= (!!(M));

        offset += XOVER::THREADS_PER_BLOCK;
    }
}

#define BIN_WIDTH 0.03125   // 1 / 32

/**
 *  This version assumes that allele_list is organized into units of WARP_SIZE (32) alleles
 *  that are ordered such that:
 *      allele_unit_idx = floor( allele_location * 32.0 )
 */
template < class StateType, class RealType >
__global__ void crossover_kernel_2(  StateType * states
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

    RealType bin_start =  (RealType) lane_id * BIN_WIDTH;

    rho *= BIN_WIDTH;   // scale rho according to number of threads per warp (bins)
    states = &states[blockIdx.x];

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
        fill_random_pool< crossover< 2 > >( states, rand_pool, bin_start, BIN_WIDTH, lane_id, (max_rounds << 5), (evt_hi - evt_lo) );
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

//            unsigned int cmask = count_preceeding_lane_events( rand_pool, _allele, lane_id, max_rounds, (evt_hi - evt_lo) );
            unsigned int cmask = count_preceeding_lane_events( rand_pool, _allele, lane_id, max_rounds );
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

/**
 *  This version assumes that allele_list is organized into units of WARP_SIZE (32) alleles
 *  that are ordered such that:
 *      allele_unit_idx = floor( allele_location * 32.0 )
 */
//template < class StateType, class RealType, unsigned int K >
template < class XOVER >
__global__ void crossover_kernel_2a( typename XOVER::state_type * states
                                    , typename XOVER::real_type * pool
                                    , typename XOVER::real_type * allele_list
                                    , unsigned int * evt_list
                                    , unsigned int * sequences
                                    , unsigned int nSequences
                                    , unsigned int nAlleles
                                    , unsigned int sequence_width
                                    , typename XOVER::poisson_type * pois_cdf ) {
    if( gridDim.x != XOVER::BLOCKS_PER_KERNEL || blockDim.x != XOVER::THREADS_PER_BLOCK ) return;

    typedef typename XOVER::state_type StateType;
    typedef typename XOVER::real_type RealType;

    typedef typename XOVER::poisson_type poisson_type;
    typedef unsigned int int_type;

    __shared__ RealType s_pois_cdf[ poisson_type::MAX_K ];
    __shared__ RealType rand_pool[ 1024 ];      // MAX_K * BINS === 32 * 32
    __shared__ int_type event_bins[ XOVER::THREADS_PER_BLOCK ];

    int_type tid = threadIdx.y * blockDim.x + threadIdx.x;
    int_type lane_id = (tid & 31);

    int_type i;

    if( tid < poisson_type::MAX_K  ) {
        s_pois_cdf[ tid ] = pois_cdf->_cdf[ tid ];
    }

    unsigned int max_k = pois_cdf->max_k;
    __syncthreads();

    RealType bin_start =  (RealType) lane_id * BIN_WIDTH;

#ifdef USE_MERSENNE_TWISTER
//    states = &states[blockIdx.x];
    __shared__ StateType local_state;

    if( tid == 0 ) {
        local_state = states[blockIdx.x];    
    }
    __syncthreads();
#else
    StateType local_state = states[ blockIdx.x * blockDim.x * blockDim.y + tid ];
#endif  // USE_MERSENNE_TWISTER

#ifdef LOG_RANDOM_EVENTS
    pool += blockIdx.x * XOVER::THREADS_PER_BLOCK;
    evt_list += blockIdx.x * XOVER::THREADS_PER_BLOCK;
#endif  // LOG_RANDOM_EVENTS

    int_type seq_idx = blockIdx.x;
    while( seq_idx < nSequences ) {

        unsigned int * seq = sequences + (seq_idx * sequence_width);

        // use all threads to generate random numbers
        RealType x = curand_uniform( &local_state );
    
        event_bins[ tid ] = _find_poisson_maxk32( s_pois_cdf, x, max_k );
        __syncthreads();

#ifdef LOG_RANDOM_EVENTS
        evt_list[ tid ] = event_bins[ tid ];
        evt_list += XOVER::BLOCKS_PER_KERNEL * XOVER::THREADS_PER_BLOCK;
#endif  // LOG_RANDOM_EVENTS

        int_type evt_lo, evt_hi;
        int_type max_rounds = warp_event_range( event_bins[lane_id], lane_id, evt_lo, evt_hi);
        __syncthreads();

        // within a warp each thread (bin) has [evt_lo, evt_hi) events

        // over generate random numbers per bin
        // avoids branch divergence
#ifdef USE_MERSENNE_TWISTER
        fill_random_pool< crossover< 2 > >( states, rand_pool, bin_start, (RealType) BIN_WIDTH, lane_id, (max_rounds << 5), (evt_hi - evt_lo) );
#else
        fill_random_pool< crossover< 2 > >( &local_state, rand_pool, bin_start, (RealType) BIN_WIDTH, lane_id, (max_rounds << 5), (evt_hi - evt_lo) );
#endif  // USE_MERSENNE_TWISTER
        __syncthreads();

#ifdef LOG_RANDOM_EVENTS
        i = tid;
        while( i < 1024 ) {
            pool[ i ] = rand_pool[ i ];
            i += THREAD_NUM;
        }

        pool += XOVER::BLOCKS_PER_KERNEL * 1024;
        __syncthreads();
#endif  // LOG_RANDOM_EVENTS

        i = tid;
        while( i < nAlleles ) {
            RealType _allele = allele_list[ i ];

            unsigned int cmask = count_preceeding_lane_events( rand_pool, _allele, lane_id, max_rounds );
            __syncthreads();

            cmask = (((cmask + evt_lo) & 1) * (1 << lane_id));  // translate event count to lane relative crossover mask

///////////////////////////////////////////////////////////
            persist_mask( cmask, (tid >> 5), lane_id, seq );
            __syncthreads();
///////////////////////////////////////////////////////////

            i += XOVER::THREADS_PER_BLOCK;
            seq += (XOVER::THREADS_PER_BLOCK >> 5);
        }
        __syncthreads();

        seq_idx += XOVER::BLOCKS_PER_KERNEL;
    }
#ifdef USE_MERSENNE_TWISTER
    if( tid == 0 ) {
        states[ blockIdx.x ] = local_state;
    }
#else
    states[ blockIdx.x * blockDim.x * blockDim.y + tid ] = local_state;
#endif  // USE_MERSENNE_TWISTER
}
crossover< 2 >::crossover( ) :
    dStates( NULL )
#ifdef USE_MERSENNE_TWISTER
    , dParams( NULL)
#endif  // USE_MERSENNE_TWISTER
{
    initialize();
}

void crossover< 2 >::initialize( ) {
    // defines from curand_mtgp32.h in curand library (as of CUDA 6.5)
    // THREAD_NUM == MTGPDC_FLOOR_2P == 256
    // BLOCK_NUM_MAX == CURAND_NUM_MTGP32_PARAMS == 200

    m_seed = clotho::utility::clock_type::now().time_since_epoch().count();

    assert( cudaMalloc( (void **) &dPoisCDF, sizeof( poisson_type ) ) == cudaSuccess);

#ifdef USE_MERSENNE_TWISTER
    assert( cudaMalloc( (void ** ) &dStates, BLOCKS_PER_KERNEL * sizeof( state_type ) ) == cudaSuccess );
    assert( cudaMalloc( (void ** ) &dParams, sizeof( state_param_type ) ) == cudaSuccess );

    assert( curandMakeMTGP32Constants( MTGPDC_PARAM_TABLE, dParams ) == CURAND_STATUS_SUCCESS );

    assert( curandMakeMTGP32KernelState( dStates, MTGPDC_PARAM_TABLE, dParams, BLOCKS_PER_KERNEL, m_seed) == CURAND_STATUS_SUCCESS );
#else

    assert( cudaMalloc( (void ** ) &dStates, BLOCKS_PER_KERNEL * THREADS_PER_BLOCK * sizeof( state_type ) ) == cudaSuccess );
    setup_state_kernel<<< BLOCKS_PER_KERNEL, THREADS_PER_BLOCK >>>( dStates, m_seed );

#endif  // USE_MERSENNE_TWISTER
}

void crossover< 2 >::adjust_alleles( allele_type * allele_list, size_type N ) {
    remap_alleles<<< 1, THREADS_PER_BLOCK >>>( allele_list, N );
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

    rho /= 32.0;

    make_poisson_cdf_maxk32<<< 1, 32 >>>( dPoisCDF, rho );

    crossover_kernel_2a< crossover< 2 > ><<< BLOCKS_PER_KERNEL, THREADS_PER_BLOCK >>>( dStates, rand_pool, allele_list, event_list, sequences, nSequences, nAlleles, sequence_width, dPoisCDF );
}

void crossover< 2 >::get_state( boost::property_tree::ptree & n ) {
    n.put( "crossover.version", 2 );
    n.put( "crossover.curand.state_type", clotho::cuda::curand_helper< typename crossover< 2 >::state_type >::StateName );
    n.put( "crossover.curand.seed", m_seed );

#ifdef USE_TEXTURE_MEMORY_FOR_ALLELE
    n.put( "crossover.use_texture_memory_for_allele", true );
#else
    n.put( "crossover.use_texture_memory_for_allele", false );
#endif  // USE_TEXTURE_MEMORY_FOR_ALLELE
}

crossover< 2 >::~crossover() {
    cudaFree( dStates );
    cudaFree( dPoisCDF );

#ifdef USE_MERSENNE_TWISTER
    cudaFree( dParams );
#endif  // USE_MERSENNE_TWISTER
}

#endif  // CROSSOVER_MATRIX_2_CUH_
