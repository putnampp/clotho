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
#ifndef RECOMBINE_PARENTS_HPP_
#define RECOMBINE_PARENTS_HPP_

#include <cuda.h>
#include "clotho/utility/algorithm_version.hpp"

template < class SequenceSpaceType, class VectorSpaceType, unsigned char V >
__global__ void recombine_parents_kernel( SequenceSpaceType * parents
                                        , VectorSpaceType * ids
                                        , SequenceSpaceType * offspring
                                        , clotho::utility::algo_version< V > * v);

#include "clotho/cuda/data_spaces/sequence_space/device_sequence_space.hpp"
#include "clotho/cuda/data_spaces/basic_data_space.hpp"

template < class IntType >
__global__ void recombine_parents_kernel( device_sequence_space< IntType > * parents
                                        , basic_data_space< IntType > * ids
                                        , device_sequence_space< IntType > * offspring 
                                        , clotho::utility::algo_version< 1 > * v) {

    typedef device_sequence_space< IntType > sequence_space_type;
    typedef typename sequence_space_type::int_type  sequence_int_type;

    typedef basic_data_space< IntType >             event_space_type;
    typedef typename event_space_type::value_type   event_int_type;

    unsigned int tpb = blockDim.x * blockDim.y;
    unsigned int bpg = gridDim.x * gridDim.y;

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;

    unsigned int parent_width = parents->seq_width;
    unsigned int off_width = offspring->seq_width;
    unsigned int off_count = offspring->seq_count;

    sequence_int_type   * par_seqs = parents->sequences;
    sequence_int_type   * off_seqs = offspring->sequences;

    event_int_type      * evts = ids->data;
    
    while(bid < off_count ) {

        unsigned int pidx = evts[bid];

        unsigned int p0_idx = pidx * parent_width + tid;
        unsigned int off_idx = bid * off_width + tid;
        unsigned int j = tid;

        int _width = ((pidx & 1) ? -parent_width : parent_width);
        while( j < parent_width ) {
            sequence_int_type p0 = par_seqs[ p0_idx ];
            sequence_int_type p1 = par_seqs[ p0_idx + _width ] ;

            sequence_int_type cross = off_seqs[ off_idx ];

            sequence_int_type off = ((p0 & ~cross) | (p1 & cross));

            off_seqs[off_idx ] = off;

            j += tpb;
            p0_idx += tpb;
            off_idx += tpb;
        }
        __syncthreads();

        // clear the tail of a offspring sequence
        while( j < off_width ) {
            off_seqs[ off_idx ] = 0;
            j += tpb;
            off_idx += tpb;
        }
        __syncthreads();

        bid += bpg;
    }
}

template < class IntType >
__global__ void recombine_parents_kernel( device_sequence_space< IntType > * parents
                                        , basic_data_space< IntType > * ids
                                        , device_sequence_space< IntType > * offspring 
                                        , clotho::utility::algo_version< 2 > * v) {

    typedef device_sequence_space< IntType > sequence_space_type;
    typedef typename sequence_space_type::int_type  sequence_int_type;

    typedef basic_data_space< IntType >             event_space_type;
    typedef typename event_space_type::value_type   event_int_type;

    unsigned int tpb = blockDim.x * blockDim.y;
    assert( (tpb & 31) == 0 );  // all warps are full

    unsigned int wpb = (tpb >> 5);
    unsigned int bpg = gridDim.x * gridDim.y;
    unsigned int spg = wpb * bpg;

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;

    unsigned int lane_id = (tid & 31);
    unsigned int warp_id = (tid >> 5 );

    unsigned int parent_width = parents->seq_width;
    unsigned int off_width = offspring->seq_width;
    unsigned int off_count = offspring->seq_count;

    sequence_int_type   * par_seqs = parents->sequences;
    sequence_int_type   * off_seqs = offspring->sequences;

    event_int_type      * evts = ids->data;

    unsigned int max_off_count = off_count / wpb;
    max_off_count += ((off_count % wpb) ? 1 : 0);
    max_off_count *= wpb;

    unsigned int off_id = (bid * wpb) + warp_id;
    
    while(off_id < max_off_count ) {    // allows blocks to terminate 'early'
        unsigned int pidx = 0;
        if( off_id < off_count ) {  // true for all threads in warp
            pidx = evts[off_id];
        }
        __syncthreads();

        unsigned int p_end = (pidx + 1) * parent_width;
        unsigned int p_start = p_end - ((off_id < off_count) ? parent_width : 0);
        p_start += lane_id;

        unsigned int end = (off_id + 1) * off_width;
        unsigned int start = end - ((off_id < off_count) ? off_width : 0);
        start += lane_id;

        int _width = ((pidx & 1) ? -parent_width : parent_width);
        while( p_start < p_end ) {
            sequence_int_type p0 = par_seqs[ p_start ];
            sequence_int_type p1 = par_seqs[ p_start + _width ] ;

            sequence_int_type cross = off_seqs[ start ];

            sequence_int_type off = ((p0 & ~cross) | (p1 & cross));

            off_seqs[ start ] = off;

            p_start += 32;
            start += 32;
        }
        __syncthreads();

        // clear the tail of a offspring sequence
        while( start < end ) {
            off_seqs[ start ] = 0;
            start += 32;
        }
        __syncthreads();

        off_id += spg;
    }
}

#endif  // RECOMBINE_PARENTS_HPP_
