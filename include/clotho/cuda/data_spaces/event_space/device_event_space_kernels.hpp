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

template < class IntType >
__global__ void _delete_space( device_event_space< IntType > * space ) {
    typename device_event_space< IntType >::int_type * evts = space->event_count;

    if( evts ) {
        delete evts;
    }
}

template < class IntType >
__global__ void _resize_space( device_event_space< IntType > * space, unsigned int N ) {
    typedef device_event_space< IntType > space_type;

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
}

#endif  // DEVICE_EVENT_SPACE_KERNELS_HPP_
