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
#ifndef SELECT_AND_CROSSOVER_UNIT_ORDERED_IMPL_HPP_
#define SELECT_AND_CROSSOVER_UNIT_ORDERED_IMPL_HPP_

#include "clotho/cuda/data_spaces/sequence_space/device_sequence_space.hpp"
#include "clotho/cuda/data_spaces/event_space/device_event_space.hpp"
#include "clotho/cuda/data_spaces/free_space/device_free_space.hpp"

#include "clotho/cuda/distributions/poisson_distribution.hpp"
#include "clotho/cuda/data_spaces/tags/unit_ordered_tag.hpp"
#include "clotho/cuda/data_spaces/tags/no_order_tag.hpp"

template < class StateType
           , class AlleleSpaceType
           , class SequenceIntType
           , class RealType
           , class IntType
        >
__global__ void select_and_crossover_kernel( StateType * states
                                            , device_sequence_space< SequenceIntType > * parent_seqs
                                            , device_event_space< IntType, no_order_tag >  * parent_ids
                                            , AlleleSpaceType * alleles
                                            , device_free_space< SequenceIntType, unit_ordered_tag > * free_space
                                            , poisson_cdf< RealType, 32 > * pois
                                            , device_sequence_space< SequenceIntType > * offspring_seqs ) {
    typedef StateType                                       state_type;
    typedef poisson_cdf< RealType, 32 >                     poisson_type;
    typedef device_sequence_space< SequenceIntType >        sequence_space_type;
    typedef typename sequence_space_type::int_type          sequence_int_type;

    typedef device_event_space< IntType, no_order_tag >    selection_type;
    typedef typename selection_type::int_type               selection_int_type;

    typedef AlleleSpaceType                                 allele_space_type;
    typedef typename allele_space_type::real_type           real_type;

    unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;
    unsigned int offspring_count = offspring_seqs->seq_count;

    if( bid >= offspring_count ) return;

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int warp_id = (tid >> 5 );
    unsigned int lane_id = (tid & 31 );

    unsigned int tpb = (blockDim.x * blockDim.y);
    unsigned int wpb = (tpb >> 5);
    unsigned int bpg = (gridDim.x * gridDim.y);
    unsigned int par_count = parent_seqs->seq_count;

    if( par_count == 0 ) {
        if( bid == 0 && tid == 0 ) {
            unsigned int off_cap = offspring_seqs->capacity;

            memset( offspring_seqs->sequences, 0, off_cap * sizeof( sequence_int_type) );
        }

        return;
    }

    state_type local_state = states[ bid * tpb + tid ];

    __shared__ unsigned int s_event_hash[ 32 ];
    __shared__ real_type    s_rand_pool[ 1024 ];
    __shared__ real_type    s_pois_cdf[ poisson_type::MAX_K ];

    unsigned int max_k = pois->max_k;
    if( tid < poisson_type::MAX_K ) {
        s_pois_cdf[ tid ] = pois->_cdf[ tid ];
    }
    __syncthreads();

    while( bid < offspring_count ) {
        real_type r = curand_uniform( &local_state );

        unsigned int e = _find_poisson_maxk32( s_pois_cdf, r, max_k );
        __syncthreads();

        if( tid < 32 ) {
            s_event_hash[ tid ] = e;
        }
        __syncthreads();

        e = s_event_hash[ lane_id ];
        __syncthreads();
        
        unsigned int psum = e;
        for( unsigned int i = 1; i < 32; i <<= 1 ) {
            unsigned int _e = __shfl_up( e, i );
            e = ((_e > e) ? _e : e );
        }

        unsigned int max = __shfl( e, 31 );

        for( unsigned int i = 1; i < 32; i <<= 1 ) {
            unsigned int p = __shfl_up( psum, i );
            psum += (( lane_id >= i ) * p );
        }

        unsigned int base_count = __shfl_up( psum, 1 );
        if( lane_id == 0 ) base_count = 0;
        __syncthreads();

        for( unsigned int i = warp_id; i < max; i += wpb ) {    // warp divergence
            r = curand_uniform( &local_state );

            s_rand_pool[ i * 32 + lane_id ] = (((real_type) lane_id + r) / (real_type) 32);
        }
        __syncthreads();

        // set up parent and offspring indexes
        selection_int_type id = par_ids[ bid ]; // all threads in block share same parent
        id <<= 1;

        unsigned int p0 = id * par_width + warp_id;
        unsigned int o_idx = bid * off_width + warp_id;
        unsigned int a_id = tid;

        while( a_id < par_allele_cap ) {    // true for all threads
            sequence_int_type p = par_seqs[ p0 ];
            sequence_int_type q = par_seqs[ p0 + par_width ];
            real_type a = all_loc[ a_id ];  // every thread reads a location just in case

            e = (psum - base_count);
            unsigned int x = base_count;
            for( unsigned int i = 0; i < max; ++i ) {   // should not diverge
                real_type y = s_rand_pool[ i * 32 + lane_id ];
                x += ((i < e) * ( y < a));
            }
            __syncthreads();

            x = ((x & 1) << lane_id);

            for( unsigned int i = 1; i < 32; i <<= 1 ) {
                unsigned int _x = __shfl_up( x, i );
                x |= ((lane_id >= i ) * _x );
            }

            if( lane_id == 31 ) {
                off_seqs[ o_idx ] = ((p & ~x) | (q & x));
            }
            __syncthreads();

            p0 += wpb;
            o_idx += wpb;
            a_id += tpb;
        }
        bid += bpg;
    }

    bid = blockIdx.y * gridDim.x + blockIdx.x;
    local_state[ bid * tpb + tid ];
}

#endif  // SELECT_AND_CROSSOVER_UNIT_ORDERED_IMPL_HPP_
