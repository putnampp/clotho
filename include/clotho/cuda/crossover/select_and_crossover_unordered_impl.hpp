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
#ifndef SELECT_AND_CROSSOVER_UNORDERED_IMPL_HPP_
#define SELECT_AND_CROSSOVER_UNORDERED_IMPL_HPP_

#include "clotho/cuda/data_spaces/device_sequence_space.hpp"
#include "clotho/cuda/data_spaces/device_event_space.hpp"
#include "clotho/cuda/data_spaces/device_free_space.hpp"

#include "clotho/cuda/data_spaces/tags/unordered_tag.hpp"
#include "clotho/cuda/data_sapces/tags/no_order_tag.hpp"

template < class StateType
           , class AlleleSpaceType
           , class SequenceIntType
           , class IntType1
           , class IntType2 >
__global__ void select_and_crossover_kernel( StateType * states
                                            , device_sequence_space< SequenceIntType > * parent_seqs
                                            , device_event_space< IntType1, no_order_tag >  * parent_ids
                                            , AlleleSpaceType * alleles
                                            , device_free_space< SequenceIntType, unordered_tag > * free_space
                                            , device_event_space< IntType2, unorderd_tag > * xover_events
                                            , device_sequence_space< SequenceIntType > * offspring_seqs ) {

    typedef StateType state_type;
    typedef device_sequence_space< SequenceIntType >        sequence_space_type;
    typedef typename sequence_space_type::int_type          sequence_int_type;

    assert( sizeof( sequence_int_type ) * 8 == 32 );

    typedef device_event_space< IntType1, no_order_tag >    selection_type;
    typedef typename selection_type::int_type               selection_int_type;

    typedef AlleleSpaceType                                 allele_space_type;
    typedef typename allele_space_type::real_type           real_type;

    unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;
    unsigned int off_count = offspring_seqs->seq_count;

    if( bid >= off_count ) return;

    assert( off_count == offspring_seqs->seq_count );

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int warp_id = (tid >> 5 );
    unsigned int lane_id = (tid & 31 );

    unsigned int tpb = (blockDim.x * blockDim.y);
    unsigned int wpb = (tpb >> 5);
    unsigned int bpg = (gridDim.x * gridDim.y);

    assert( tpb & 31 == 0 );    // assert all warps are full

    state_type local_state = states[ bid * tpb + tid ];

    sequence_int_type   * par_seqs = parent_seqs->sequences;
    sequence_int_type   * off_seqs = offspring_seqs->sequences;

    selection_int_type  * par_ids = parent_ids->event_count;

    unsigned int par_count = parent_seqs->seq_count;
    unsigned int par_width = parent_seqs->seq_width;

    unsigned int off_width = offspring_seqs->seq_width;

    unsigned int par_allele_cap = alleles->capacity;
    unsigned int par_allele_size = alleles->size;
    real_type * all_loc = alleles->locations;

    assert( tpb == allele_space_type::ALIGNMENT_SIZE );
    assert( par_allele_cap % tpb == 0 ); // assert_1: that parent sequence width is a multiple of threads

    __shared__ real_type    s_pois_cdf[ poisson_type::MAX_K ];
    __shared__ real_type    s_rand_pool[ allele_space_type::ALIGNMENT_SIZE ];
    __shared__ unsigned int s_event_hash[ allele_space_type::ALIGNMENT_SIZE + 1 ];

    unsigned int max_k = pois->max_k;
    if( tid < poisson_type::MAX_K ) {
        s_pois_cdf[ tid ] = pois->_cdf[ tid ];
    }
    __syncthreads();

    while( bid < off_count ) {  // true for all threads in block
        // generate crossover events
        //
        
        s_event_hash[ tid ] = 0;    // clear hash

        real_type y = curand_uniform( &local_state );
        s_rand_pool[ tid ] = curand_uniform( &local_state );

        unsigned int e = _find_poisson_maxk32( s_pois_cdf, y, max_k );
        __syncthreads();

        for( unsigned int i = 1; i < 32; i <<= 1 ) {
            unsigned int tmp = __shfl_up( e, i );
            e += ((lane_id >= i) * tmp);
        }

        if( wpb > 1 ) { // true for all threads in block (no diverge)
            if( lane_id == 31 ) {
                s_event_hash[ 32 + warp_id ] = e;
            }
            __syncthreads();

            unsigned int psum = s_event_hash[ 32 + lane_id ];
            __syncthreads();
            for( unsigned int i = 1; i < 32; i <<= 1 ) {
                unsigned int tmp = __shfl_up( psum, i );
                psum += ((lane_id >= i) * tmp );
            }

            unsigned int t = 0;
            if( warp_id > 0 ) { // no diverge within warp
                t = __shfl( psum, warp_id - 1 );
            }
            __syncthreads();

            e += t;
        }

        s_event_hash[ tid + 1 ] = e;
        __syncthreads();

        unsigned int s = s_event_hash[ tid ];
        __syncthreads();

        // adjust rand_pool such that rand_pool[ [0, s_event_hash[ allele_space_type::ALIGNMENT_SIZE ] ) ] are ordered
        real_type accum = 0.;
        while( s < e ) {    // threads diverge
            y = rand_pool[ s ];

            accum += ( log(y) / (real_type)( e - s) );
            rand_pool[ s++ ] = ( (((real_type)tid) + (1.0 - exp(accum) )) / ((real_type)allele_space_type::ALIGNMENT_SIZE );
        }
        __syncthreads();

        // set up parent and offspring indexes
        selection_int_type id = par_ids[ bid ]; // all threads in block share same parent

        id <<= 1;
        assert( (id + 1) < par_count );

        unsigned int p0 = id * par_width + warp_id;
        unsigned int p1 = p0 + par_width;

        unsigned int o_idx = bid * off_width + warp_id;
        unsigned int a_id = tid;

        while( a_id < par_allele_cap ) {    // true for all threads (no divergence) otherwise assert_1 would have failed
            sequence_int_type p = par_seqs[ p0 ];
            sequence_int_type q = par_seqs[ p1 ];

            sequence_int_type x = ((a_id < par_allele_size) * ((p ^ q) & (1 << lane_id )));
            if( x ) {   // thread divergence
                real_type a = all_loc[ a_id ];

                // count prior crossover events
                unsigned int tmp = (unsigned int)( a * ((real_type)allele_space_type::ALIGNMENT_SIZE));

                s = s_event_hash[ tmp++ ];
                e = s_event_hash[ tmp ];

                tmp = s;
                while( s < e ) {
                    y = rand_pool[ s++ ];
                    tmp += ( a > y );
                }

                x *= (tmp & 1); // if odd number of events then crossover
            }
            __syncthreads();

            // reduce crossover mask x
            for( unsigned int i = 1; i < 32; i <<= 1 ) {
                sequence_int_type tmp = __shfl_up( x, i );
                x |= ( lane_id >= i )* tmp;
            }

            sequence_int_type o = ((p & ~x) | (q & x));
            if( lane_id == 31 ) {
                off_seqs[ o_idx ] = o;
            }
            __syncthreads();
            p0 += wpb;
            p1 += wpb;

            o_idx += wpb;

            a_id += tpb;
        }
        bid += bpg;
    }

    states[ bid * tpb + tid ] = local_state;
}

#endif  // SELECT_AND_CROSSOVER_UNORDERED_IMPL_HPP_
