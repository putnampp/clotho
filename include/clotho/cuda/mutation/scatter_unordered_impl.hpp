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
#ifndef SCATTER_UNORDERED_IMPL_HPP_
#define SCATTER_UNORDERED_IMPL_HPP_

#include "clotho/cuda/data_spaces/free_space/device_free_space.hpp"
#include "clotho/cuda/data_spaces/sequence_space/device_sequence_space.hpp"
#include "clotho/cuda/data_spaces/event_space/device_event_space.hpp"
#include "clotho/cuda/data_spaces/tags/unordered_tag.hpp"
#include "clotho/cuda/index_converter.hpp"

/*
template < class StateType, class RealType, class IntType >
__global__ void _scatter_mutation( StateType * states
                                    , device_allele_space< RealType, IntType, unordered_tag > * alleles
                                    , device_event_space< IntType, unordered_tag >  * events
                                    , device_sequence_space< IntType >              * sequences ) {

    typedef device_event_space< IntType, unordered_tag >            event_space_type;
    typedef typename event_space_type::int_type                     event_int_type;

    typedef device_sequence_space< IntType >    sequence_type;
    typedef typename sequence_type::int_type    sequence_int_type;

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    event_int_type *   evt_ptr = events->event_count;

    sequence_int_type    * seq_ptr = sequences->sequences;
    unsigned int N = sequences->seq_count;
    unsigned int width = sequences->seq_width;

    unsigned int padded_N = N % (blockDim.x * blockDim.y);
    padded_N = N - padded_N;    // remove tail

    unsigned int i = tid;

    fixed_width_converter< sequence_type::OBJECTS_PER_INT > converter;

    unsigned int e_min = 0;
    while( i < padded_N ) { // all threads are active
        event_int_type e_max = evt_ptr[ i ];

        event_int_type t = __shfl_up( e_max, 1 );
        e_min = ((tid == 0 ) ? e_min : t );

        // begin divergent code
        // each active thread updates global memory
        // threads access memory in different strides
        // stride offset are not uniform
        //
        // Is the divergence cost acceptable?
        while( e_min < e_max ) {
            unsigned int idx = free_map[e_min++];   // lookup free index

            unsigned int block_idx = converter.major_offset( idx );
            idx = converter.minor_offset( idx );

            sequence_int_type b = seq_ptr[ i * width + block_idx ];
            b |= ( 1 << idx );
            seq_ptr[ i * width + block_idx ] = b;
        }
        __syncthreads();
        // end divergent code

        i += (blockDim.x * blockDim.y);

        e_min = __shfl( e_max, 31 );
    }

    if( i < N ) {   // some threads are active

    }
}*/

template < class IntType >
__global__ void _scatter_mutation_single_thread( device_free_space< IntType, unordered_tag > * fspace 
                                    , device_event_space< IntType, unordered_tag >  * events
                                    , device_sequence_space< IntType >              * sequences
                                    , unsigned int offset ) {

    typedef device_free_space< IntType, unordered_tag >             free_space_type;

    typedef device_event_space< IntType, unordered_tag >            event_space_type;
    typedef typename event_space_type::int_type                     event_int_type;

    typedef device_sequence_space< IntType >    sequence_type;
    typedef typename sequence_type::int_type    sequence_int_type;

//    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x + offset;

    event_int_type  e_max = events->event_count[ bid ];
    event_int_type  e_min = ((bid == 0 ) ? 0 : events->event_count[ bid - 1]);

    __syncthreads();
    
    if( e_min >= e_max ) return;

    fixed_width_converter< sequence_type::OBJECTS_PER_INT > converter;

    free_space_type loc_free = *fspace;

    sequence_int_type seq_width = sequences->seq_width;
    sequence_int_type * seq = sequences->sequences;
    seq += bid * seq_width;

    unsigned int cur_block = 0;
    sequence_int_type b = seq[cur_block]; // read current bit_block

    while( e_min < e_max ) {
        unsigned int fidx = loc_free.free_map[ e_min++ ];

        unsigned int block_idx = converter.major_offset( fidx );
        unsigned int bit_idx = converter.minor_offset( fidx );

        if( block_idx != cur_block ) {
            seq[cur_block] = b;   // persist current bit block to

            cur_block = block_idx;             // update current bit block

            b = seq[cur_block];   // read current bit block
        }

        b |= ((sequence_int_type) 1 << bit_idx );
    }


    seq[cur_block] = b;   // persist current bit block
}
                                
#endif  // SCATTER_UNORDERED_IMPL_HPP_
