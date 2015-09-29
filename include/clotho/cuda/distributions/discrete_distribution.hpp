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
#ifndef CLOTHO_CUDA_DISCRETE_DISTRIBUTION_HPP_
#define CLOTHO_CUDA_DISCRETE_DISTRIBUTION_HPP_

#include <cuda.h>
#include "clotho/cuda/distributions/discrete_table.hpp"
#include "clotho/cuda/popcount_kernel.h"

#include "clotho/cuda/device_state_object.hpp"

/****************************************************************
 *
 * The following algorithms are based upon the discrete_distribution
 * provided by the Boost Random library (v 1.56.0)
 *
 ***************************************************************/

/**
 *
 * Compute the average of the weights
 * Normalize each weight relative to the average
 *
 */
template < class RealType, class IntType >
__device__ void init_weights( RealType * in, RealType * out, IntType * indices, unsigned int N, unsigned int id ) {
    typedef RealType real_type;

    // set block relative thread id
    unsigned int tid = id;
    unsigned int warp_id = (tid >> 5);  // assumes 32 threads/warp
    unsigned int lane_id = (tid & 31);  // assumes 32 threads/warp

    unsigned int tcount = blockDim.x * blockDim.y;

    real_type _s = 0.;

    while( tid < N ) {
        real_type v = in[tid];
        _s += v;
        tid += tcount;
    }
    __syncthreads();

    // reset thread id
    tid = id;

    for( unsigned int i = 1; i < 32; i <<= 1 ) {
        real_type s = __shfl_up( _s, i );
        _s += ((real_type)( lane_id >= i )) * s;
    }

    if( tcount > 32 ) { // true/false for all threads in block
        __shared__ real_type buffer[ 32 ];  // max of 32 warps/block

        if( warp_id == 0 ) {    // use the threads in warp_0 to clear buffer 
            buffer[lane_id] = 0.;
        }
        __syncthreads();

        buffer[ warp_id ] = _s;
        __syncthreads();

        _s = buffer[lane_id];

        for( unsigned int i = 1; i < 32; i <<= 1 ) {
            real_type s = __shfl_up( _s, i );
            _s += ((real_type)( lane_id >= i )) * s;
        }
    }
 
    // every warp computes the total summation into thread 31 of warp   
    real_type ave = __shfl( _s, 31 );
    ave /= ((real_type) N);

    while( tid < N ) {
        real_type v = in[tid];
        v /= ave;

        out[tid] = v;
        indices[tid] = tid;
        tid += tcount;
    }
}

template < class IntType, class RealType >
__global__ void init_discrete_distribution( basic_data_space< RealType > * weights
                                            , discrete_table< IntType, RealType > * tbl
                                            , discrete_table< IntType, RealType > * above_buffer
                                            , discrete_table< IntType, RealType > * below_buffer ) {

    typedef discrete_table< IntType, RealType > table_type;
    typedef typename table_type::real_type      real_type;
    typedef typename table_type::int_type       int_type;

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int N = weights->size;
    if( tid == 0 ) {
        _resize_space_impl( tbl, N );
        _resize_space_impl( above_buffer, N );
        _resize_space_impl( below_buffer, N );
    }
    __syncthreads();

    init_weights( weights->data, tbl->threshold, tbl->alternative, N, tid);
}

template < class IntType, class RealType >
__global__ void partition_discrete_distribution( discrete_table< IntType, RealType > * tbl
                                                , discrete_table< IntType, RealType > * above_buffer
                                                , discrete_table< IntType, RealType > * below_buffer ) {

    typedef discrete_table< IntType, RealType > table_type;
    typedef typename table_type::real_type      real_type;
    typedef typename table_type::int_type       int_type;
    
    unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;
    if( bid >= 2 ) return;

    real_type   * thresh_base = tbl->threshold;
    int_type    * alt_base = tbl->alternative;

    real_type   * thresh;
    int_type    * alt;

    if( bid == 0 ) {
        thresh = above_buffer->threshold;
        alt = above_buffer->alternative;
    } else {
        thresh = below_buffer->threshold;
        alt = below_buffer->alternative;
    }

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int lane_id = (tid & 31);
    unsigned int N = tbl->size;

    popcountGPU< unsigned int > pc;
    unsigned int offset = 0;
    unsigned int tail = (N % (blockDim.x * blockDim.y));

    unsigned int paddedN = (N - tail) + (blockDim.x * blockDim.y);

    while( tid < paddedN ) {
        real_type t = 0;
        int_type a = 0;
        
        if( tid < N ) {
            t = thresh_base[ tid ];
            a = alt_base[ tid ];
        } else {
            t = ((bid == 0) ? 0. : 1.);
            a = 0;
        }
        __syncthreads();

        unsigned int mask = (1 << lane_id);
        if( bid == 0 ) {
            mask *= (!( t < 1.0 ));
        } else {
            mask *= ( t < 1.0 );
        }

        for( unsigned int i = 1; i < 32; i <<= 1 ) {
            unsigned int m = __shfl_up( mask, i );
            mask |= (( lane_id >= i ) * m);
        }

        unsigned int m = __shfl( mask, 31 );    // warp mask

        mask = (1 << lane_id );
        unsigned int idx = pc.evalGPU( m & (mask - 1));
        if( m & mask && (tid < N) ) {
            // move elements
            //
            idx += offset;
            thresh[idx] = t;
            alt[idx] = a;
        }
        __syncthreads();

        offset += pc.evalGPU( m );
        tid += (blockDim.x * blockDim.y);
    }

    if( threadIdx.y * blockDim.x + threadIdx.x == 0 ) {
        if( bid == 0 ) {
            above_buffer->size = offset;
        } else {
            below_buffer->size = offset;
        }
    }
}


/**
 *
 * Uncertain of how to best parallelize this algorithm
 *
 */
/*
template < class IntType, class RealType >
__global__ void finalize_discrete_distribution( discrete_table< IntType, RealType > * tbl
                                                , discrete_table< IntType, RealType > * above_buffer
                                                , discrete_table< IntType, RealType > * below_buffer ) {

    typedef discrete_table< IntType, RealType > table_type;
    typedef typename table_type::real_type      real_type;
    typedef typename table_type::int_type       int_type;

    unsigned int N = tbl->size;
    unsigned int A = above_buffer->size;
    unsigned int B = below_buffer->size;

    real_type * thresh_base = tbl->threshold;
    real_type * a_thresh = above_buffer->threshold;
    real_type * b_thresh = below_buffer->threshold;

    int_type * alt_base = tbl->alternative;
    int_type * a_alt = above_buffer->alternative;
    int_type * b_alt = below_buffer->alternative;

    assert( A + B == N );

    unsigned int a = 0, b = 0;

    while( b < B && a < A ) {
        int_type    b_idx = b_alt[ b ];
        int_type    a_idx = a_alt[a];

        real_type b_val = b_thresh[ b ];
        real_type a_val = a_thresh[ a ];

        thresh_base[ b_idx ] = b_val;
        alt_base[ b_idx ] = a_idx;

        a_val -= (1.0 - b_val);

        if( a_val < 1.0 ) {
            b_thresh[b] = a_val;
            b_alt[b] = a_idx;
            ++a;
        } else {
            ++b;
        }
    }

    while( b < B ) {
        int_type idx = b_alt[ b++ ];
        thresh_base[ idx ] = 1.0;
    }

    while( a < A ) {
        int_type idx = a_alt[ a++ ];

        thresh_base[ idx ] = 1.0;
    }
}*/

/**
 *
 * SINGLE BLOCK SINGLE WARP parallelization. 
 *
 */
template < class IntType, class RealType >
__global__ void finalize_discrete_distribution2( discrete_table< IntType, RealType > * tbl
                                                , discrete_table< IntType, RealType > * above_buffer
                                                , discrete_table< IntType, RealType > * below_buffer ) {

    typedef discrete_table< IntType, RealType > table_type;
    typedef typename table_type::real_type      real_type;
    typedef typename table_type::int_type       int_type;

    unsigned int N = tbl->size;
    unsigned int A = above_buffer->size;
    unsigned int B = below_buffer->size;

    real_type * thresh_base = tbl->threshold;
    real_type * a_thresh = above_buffer->threshold;
    real_type * b_thresh = below_buffer->threshold;

    int_type * alt_base = tbl->alternative;
    int_type * a_alt = above_buffer->alternative;
    int_type * b_alt = below_buffer->alternative;

    assert( A + B == N );

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int lane_id = (tid & 31);

    __shared__ unsigned int a_stack[ 32 ];
    __shared__ unsigned int b_stack[ 32 ];

    a_stack[ tid ] = tid;
    b_stack[ tid ] = tid;

    popcountGPU< unsigned int > pc;

    while( true ) {

        unsigned int a = a_stack[ tid ];
        unsigned int b = b_stack[ tid ];
        __syncthreads();

        // horribly divergent code
        unsigned int m = 0, all_done = 0;
        if( a >= A ) {
            if( b < B ) {
                int_type _idx = b_alt[ b ];
                thresh_base[ _idx ] = 1.0;
            } else {
                all_done = (1 << lane_id);
            }
        } else if( b >= B ) {
            int_type _idx = a_alt[ a ];
            thresh_base[ _idx ] = 1.0;
            m = ( 1 << lane_id );
        } else {
            int_type b_idx = b_alt[ b ];
            int_type a_idx = a_alt[ a ];

            real_type   b_val = b_thresh[ b ];
            real_type   a_val = a_thresh[ a ];

            thresh_base[ b_idx ] = b_val;
            alt_base[ b_idx ] = a_idx;

            a_val -= (1.0 - b_val);

            if( a_val < 1.0 ) { // remaining complement of current above element is now below normalized average (1.0)
                // move complement to below stack
                b_thresh[ b ] = a_val;
                b_alt[ b ] = a_idx;
                m = (1 << lane_id );
            } else {
                a_thresh[ a ] = a_val;
            }
        }
        __syncthreads();
        // end of horribleness

        for( unsigned int i = 1; i < 32; i <<= 1 ) {
            unsigned int _d = __shfl_up( all_done, i );
            all_done |= ((lane_id >= i) * _d);
        }

        unsigned int _d = __shfl( all_done, 31 );
        if( (~_d) == 0 ) { // true for all threads
            break;
        }

        for( unsigned int i = 1; i < 32; i <<= 1 ) {
            unsigned int _m = __shfl_up( m, i );
            m |= (( lane_id >= i ) * _m);
        }

        unsigned int mask = __shfl( m, 31 );    // all 1's should update their a index;

        unsigned int a_updates = pc.evalGPU( mask );
        unsigned int b_updates = 32 - a_updates;

        __syncthreads();

        if( b_updates ) {   // true for all threads in block (no divergence)
            // if there b updates then shift the non-updates down
            //
            // find the largest index in the current set of indices
            unsigned int bmax = b;
            for( unsigned int i = 1; i < 32; i <<= 1 ) {
                unsigned int _m = __shfl_up( bmax, i );
                bmax = ((bmax < _m ) ? _m : bmax);
            }

            unsigned int _b = __shfl( bmax, 31 );

            unsigned int idx_offset = pc.evalGPU( (~mask) & (( 1 << lane_id) - 1) );

            if( mask & (1 << lane_id ) ) {  // all above updates are paired with below values which should be shifted down
                b_stack[ lane_id - idx_offset ] = b;
            }
            __syncthreads();

            if( lane_id >= 32 - b_updates ) {
                unsigned int tmp = _b + (b_updates - (32 - lane_id)) + 1;
                b_stack[ lane_id ] =  ((tmp < B) ? tmp : B);
            }
            __syncthreads();
        }

        if( a_updates ) {   // true for all threads in block (no divergence)

            unsigned int amax = a;
            for( unsigned int i = 1; i < 32; i <<= 1 ) {
                unsigned int _m = __shfl_up( amax, i );
                amax = (( amax < a ) ? a : amax);
            }

            unsigned int _a = __shfl( amax, 31 );
            unsigned int idx_offset = pc.evalGPU( (mask) & (( 1 << lane_id ) - 1));

            if( (~mask) & (1 << lane_id ) ) {
                a_stack[ lane_id - idx_offset ] = a;
            }
            __syncthreads();

            if( lane_id >= (32 - a_updates) ) {
                unsigned int tmp =_a + ( a_updates - (32 - lane_id ) + 1 );
                a_stack[ lane_id ] = ((tmp < A) ? tmp : A);
            }
            __syncthreads();

        }
    }
}

template < class StateType, class IntType, class RealType >
__global__ void make_random_list( StateType * states
                                  , discrete_table< IntType, RealType > * tbl
                                  , device_sequence_space< IntType > * parents
                                  , device_sequence_space< IntType > * offspring
                                  , device_event_space< IntType, no_order_tag > * events ) {
    typedef StateType state_type;

    typedef discrete_table< IntType, RealType > table_type;
    typedef typename table_type::int_type       int_type;
    typedef typename table_type::real_type      real_type;

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    state_type local_state = states[ tid ];

    unsigned int nParents = parents->seq_count / 2;
    unsigned int nOff = offspring->seq_count;

    if( tid == 0 ) {
        _resize_space_impl( events, nOff );
    }
    __syncthreads();

    unsigned int t_size = tbl->size;

    assert( t_size == nParents );   // true for all threads

    unsigned int * event_counts = events->event_count;
    unsigned int _count = blockDim.x * blockDim.y;
    if( t_size == 0 ) {
        while( tid < nOff ) {
            event_counts[ tid ] = 0;
            tid += _count;
        }
        return;   
    }

    __syncthreads();

    real_type * tlist = tbl->threshold;
    int_type *  alts = tbl->alternative;

    unsigned int nPad = nOff % _count;
    nPad = (nOff - nPad) + _count;

    while( tid < nPad ) {  
        real_type p = curand_uniform( &local_state );
        int_type  idx = curand( &local_state );
        idx %= nParents;

        real_type q = tlist[ idx ];
        int_type tmp = alts[ idx ];
        
        idx = ((p < q) ? idx : tmp);

        if( tid < nOff ) {
            event_counts[ tid ] = idx;
        }
        __syncthreads();
        tid += _count;
    }

    states[ threadIdx.y * blockDim.x + threadIdx.x ] = local_state;
}

template < class IntType, class RealType >
class DiscreteDistribution : public clotho::utility::iStateObject {
public:
    typedef discrete_table< IntType, RealType > table_type;

    typedef basic_data_space< RealType >        weight_space_type;

    DiscreteDistribution() :
        dTable( NULL )
        , dAboveBuffer( NULL )
        , dBelowBuffer( NULL )
    {
        create_space( dTable );
        create_space( dAboveBuffer );
        create_space( dBelowBuffer );
    }

    void initialize_table( weight_space_type * weights ) {
        init_discrete_distribution<<< 1, 32 >>>( weights, dTable, dAboveBuffer, dBelowBuffer );
        CHECK_LAST_KERNEL_EXEC

        partition_discrete_distribution<<< 2, 32 >>>( dTable, dAboveBuffer, dBelowBuffer );
        CHECK_LAST_KERNEL_EXEC

        finalize_discrete_distribution2<<< 1, 32 >>>( dTable, dAboveBuffer, dBelowBuffer );
        CHECK_LAST_KERNEL_EXEC
    }

    table_type *  get_device_space() {
        return dTable;
    }

    void get_state( boost::property_tree::ptree & state ) {
        boost::property_tree::ptree tbl;
        get_device_object_state( tbl, dTable );
        state.add_child( "table", tbl );

        /*
        boost::property_tree::ptree abv, bel;
        get_device_state_object( abv, dAboveBuffer );
        get_device_state_object( bel, dBelowBuffer );

        state.add_child( "above", abv );
        state.add_child( "below", bel );
        */

    }

    virtual ~DiscreteDistribution() {
        delete_space( dBelowBuffer );
        delete_space( dAboveBuffer );
        delete_space( dTable );
    }
protected:
    table_type  * dTable;
    table_type  * dAboveBuffer, * dBelowBuffer;
};

#endif  // CLOTHO_CUDA_DISCRETE_DISTRIBUTION_HPP_
