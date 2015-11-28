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
#ifndef DEVICE_EVENT_SPACE_KERNELS_HPP_
#define DEVICE_EVENT_SPACE_KERNELS_HPP_

#include "clotho/cuda/data_spaces/event_space/device_event_space_def.hpp"
#include "clotho/cuda/data_spaces/data_space_kernel_api.hpp"


template < class IntType, class OrderTag >
__global__ void _delete_space( device_event_space< IntType, OrderTag > * space ) {
//    assert( blockIdx.y * gridDim.x + blockIdx.x == 0 );
//    assert( threadIdx.y * blockDim.x + threadIdx.x == 0 );
//
//    typename device_event_space< IntType, OrderTag >::int_type * evts = space->event_count;
//    space->event_count = NULL;
//
//    if( evts ) {
//        delete evts;
//    }
}

/*
template < class IntType, class OrderTag >
__device__ void _resize_space_impl( device_event_space< IntType, OrderTag > * space, unsigned int N ) {
    typedef device_event_space< IntType, OrderTag > space_type;

    if( space->capacity < N ) {
        typename space_type::pointer tmp = space->event_count;
        if( tmp ) {
            delete tmp;
        }

        space->event_count = new typename space_type::int_type[ N ];
        space->capacity = N;
    }

    space->size = N;
    space->total = 0;
}*/

template < class IntType, class OrderTag >
__global__ void _resize_space( device_event_space< IntType, OrderTag > * space, unsigned int N ) {
//    assert(blockIdx.y * gridDim.x + blockIdx.x == 0 );
//
//    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
//    assert( tid == 0);
//
//    if( tid == 0 ) {
//        _resize_space_impl( space, N );
//    }
}

#endif  // DEVICE_EVENT_SPACE_KERNELS_HPP_
