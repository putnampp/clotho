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
#ifndef RANDOM_SAMPLE_SEQUENCE_SPACE_HPP_
#define RANDOM_SAMPLE_SEQUENCE_SPACE_HPP_

#include "clotho/cuda/sampling/random_sample_def.hpp"

#include "clotho/cuda/data_spaces/sequence_space/device_sequence_space.hpp"
#include "clotho/cuda/data_spaces/basic_data_space.hpp"

template < class StateType, class IntType, class IntType2 = unsigned int >
__global__ void random_select_parents_kernel( StateType * states
                                                , device_sequence_space< IntType > * src
                                                , unsigned int N
                                                , basic_data_space< IntType2 > * index_space ) {
    typedef StateType   state_type;

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    state_type  local_state = states[ tid ];

    IntType2 M = src->seq_count;

    if( tid == 0 ) {
        _resize_space_impl( index_space, N );
    }
    __syncthreads();

    unsigned int i = tid;
    unsigned int * ids = index_space->data;
    unsigned int tpb = blockDim.x * blockDim.y;

    while( i < N ) {
        float x = curand_uniform( &local_state );

        IntType2 pidx = (IntType2)(x * M);

        pidx = ((pidx >= M) ? 0 : pidx);    // wrap around indexes

        ids[ i ] = pidx;

        i += tpb;
    }

    states[tid] = local_state;
}
#endif  // RANDOM_SAMPLE_SEQUENCE_SPACE_HPP_
