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
#ifndef DEVICE_EVENT_POOL_SPACE_HPP_
#define DEVICE_EVENT_POOL_SPACE_HPP_

#include "clotho/cuda/data_spaces/data_space_kernel_api.hpp"

template < class EventType, class IndexOffsetType >
struct device_event_pool_space {
    typedef EventType   event_type;
    typedef IndexOffsetType index_offset_type;

    static const unsigned int MAX_EVENTS_PER_INDEX = 32;

    event_type          * event_pool;
    index_offset_type   * offsets;

    unsigned int pool_size, pool_capacity;
    unsigned int offset_size, offset_capacity;
};

template < class EventType, class IndexOffsetType >
__global__ void _resize_space( device_event_pool_space< EventType, IndexOffsetType > * space, unsigned int N ) {
    assert(blockIdx.y * gridDim.x + blockIdx.x == 0 );

    if( threadIdx.y * blockDim.x + threadIdx.x == 0 ) {
        typedef device_event_pool_space< EventType, IndexOffsetType > space_type;

        unsigned int pool_cap = space->pool_capacity;
        unsigned int new_pool_size = space_type::MAX_EVENTS_PER_OFFSET * N;
        if( pool_cap < new_pool_size ) {
            if( space->event_pool ) {
                delete [] space->event_pool;
            }
            space->event_pool = new typename space_type::event_type[ new_pool_size ];
            space->pool_capacity = new_pool_size;
        }
        space->pool_size = new_pool_size;

        if( space->offset_capacity < N + 1 ) {
            if( space->offsets ) {
                delete [] space->offsets;
            }

            space->offsets = new typename space_type::index_offset_type[ N + 1 ];
            space->offset_capacity = N + 1;
        }
        space->offset_capacity = N + 1;
    }
}


template < class EventType, class IndexOffsetType >
__global__ void _delete_space( device_event_pool_space< EventType, IndexOffsetType > * space, unsigned int N ) {
    if( blockIdx.y * gridDim.x + blockIdx.x != 0 ) return;

    if( threadIdx.y * blockDim.x + threadIdx.x != 0 ) return;

    typename device_event_pool_space< EventType, IndexOffsetType >::event_type * tmp_pool = space->event_pool;
    if( tmp_pool ) {
        delete [] tmp_pool;
    }

    space->event_pool = NULL;

    typename device_event_pool_space< EventType, IndexOffsetType >::index_offset_type * tmp_off = space->offsets;
    if( tmp_off ) {
        delete [] tmp_off;
    }

    space->offsets = NULL;
}


#endif  // DEVICE_EVENT_POOL_SPACE_HPP_
