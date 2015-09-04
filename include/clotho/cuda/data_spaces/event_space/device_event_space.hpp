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
#ifndef DEVICE_EVENT_SPACE_HPP_
#define DEVICE_EVENT_SPACE_HPP_

#include "clotho/cuda/data_spaces/event_space/device_event_space_def.hpp"
#include "clotho/cuda/data_spaces/event_space/device_event_space_kernels.hpp"

/*
template < class IntType, class OrderTag >
void create_space( device_event_space< IntType, OrderTag > *& space, unsigned int N ) {
    unsigned int n = sizeof( device_event_space< IntType, OrderTag > );
    assert( cudaMalloc( (void **) &space, n ) == cudaSuccess );
    assert( cudaMemset( space, 0, n ) == cudaSuccess );

    if( N ) {
        resize_space( space, N );
    }
}

template < class IntType, class OrderTag >
void resize_space( device_event_space< IntType, OrderTag > * space, unsigned int N ) {
    _resize_space<<< 1, 1 >>>( space, N );
}

template < class IntType, class OrderTag >
void delete_space( device_event_space< IntType, OrderTag > * space ) {
    _delete_space<<< 1, 1 >>>( space );

    cudaFree( space );
}*/

#endif  // DEVICE_EVENT_SPACE_HPP_
