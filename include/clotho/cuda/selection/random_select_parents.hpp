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
#ifndef RANDOM_SELECT_PARENTS_HPP_
#define RANDOM_SELECT_PARENTS_HPP_

#include <cuda.h>

template < class StateType, class SequenceSpaceType, class EventSpaceType >
__global__ void random_select_parents_kernel( StateType * states
                                              , SequenceSpaceType * parents
                                              , SequenceSpaceType * child
                                              , EventSpaceType * events );


#include "clotho/cuda/data_spaces/sequence_space/device_sequence_space.hpp"
#include "clotho/cuda/data_spaces/event_space/device_event_space.hpp"
#include "clotho/cuda/data_spaces/tags/no_order_tag.hpp"

template < class StateType, class IntType >
__global__ void random_select_parents_kernel( StateType * states
                                                , device_sequence_space< IntType > * parents
                                                , device_sequence_space< IntType > * child 
                                                , device_event_space< IntType, no_order_tag > * events ) {

    typedef StateType   state_type;

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    state_type  local_state = states[ tid ];

    unsigned int nParents = parents->seq_count / 2;
    unsigned int nChild = child->seq_count;

    if( tid == 0 ) {
        _resize_space_impl( events, nChild );
    }
    __syncthreads();

    unsigned int i = tid;

    unsigned int * event_counts = events->event_count;

    unsigned int offset = blockDim.x * blockDim.y;

    while( i < nChild ) {
        float x = curand_uniform( &local_state );
        unsigned int pidx = (x * nParents);

        event_counts[ i ] = (pidx << 1);

        i += offset;
    }

    states[tid] = local_state;
}

#endif  // RANDOM_SELECT_PARENTS_HPP_
