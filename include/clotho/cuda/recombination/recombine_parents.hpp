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

template < class SequenceSpaceType, class EventSpaceType >
__global__ void recombine_parents_kernel( SequenceSpaceType * parents
                                        , EventSpaceType * events
                                        , SequenceSpaceType * offspring );

#include "clotho/cuda/data_spaces/sequence_space/device_sequence_space.hpp"
#include "clotho/cuda/data_spaces/event_space/device_event_space.hpp"
#include "clotho/cuda/data_spaces/tags/no_order_tag.hpp"

template < class IntType >
__global__ void recombine_parents_kernel( device_sequence_space< IntType > * parents
                                        , device_event_space< IntType, no_order_tag > * events
                                        , device_sequence_space< IntType > * offspring ) {

    typedef device_sequence_space< IntType > sequence_space_type;
    typedef typename sequence_space_type::int_type  sequence_int_type;

    typedef device_event_space< IntType, no_order_tag > event_space_type;
    typedef typename event_space_type::int_type event_int_type;

    unsigned int thread_offset = blockDim.x * blockDim.y;
    unsigned int off_offset = gridDim.x * gridDim.y;

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;

    unsigned int parent_width = parents->seq_width;
    unsigned int off_width = offspring->seq_width;
    unsigned int off_count = offspring->seq_count;

    sequence_int_type   * par_seqs = parents->sequences;
    sequence_int_type   * off_seqs = offspring->sequences;

    event_int_type      * evts = events->event_count;
    
    unsigned int i = bid;
    while(i < off_count ) {

        unsigned int pidx = evts[i];

        unsigned int p0_idx = pidx * parent_width + tid;
        unsigned int off_idx = i * off_width + tid;
        unsigned int j = tid;
        while( j < parent_width ) {
            sequence_int_type p0 = par_seqs[ p0_idx ];
            sequence_int_type p1 = par_seqs[ p0_idx + parent_width ] ;

            sequence_int_type cross = off_seqs[ off_idx ];

            sequence_int_type off = ((p0 & ~cross) | (p1 & cross));

            off_seqs[off_idx ] = off;

            j += thread_offset;
            p0_idx += thread_offset;
            off_idx += thread_offset;
        }
        __syncthreads();

        while( j < off_width ) {
            off_seqs[ off_idx ] = 0;
            j += thread_offset;
            off_idx += thread_offset;
        }
        __syncthreads();

        i += off_offset;
    }
}

#endif  // RECOMBINE_PARENTS_HPP_
