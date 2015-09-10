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
#ifndef DEVICE_FREE_SPACE_KERNELS_HPP_
#define DEVICE_FREE_SPACE_KERNELS_HPP_

#include "clotho/cuda/data_spaces/free_space/device_free_space_def.hpp"
#include "clotho/cuda/data_spaces/data_space_kernel_api.hpp"

template < class IntType, class OrderTag >
__device__ void _resize_space_impl( device_free_space< IntType, OrderTag > * s, unsigned int N ) {
    typedef device_free_space< IntType, OrderTag > space_type;

    unsigned int free_list_size = N / space_type::OBJECTS_PER_INT;
    if( N % space_type::OBJECTS_PER_INT ) {
        ++free_list_size;
    }

    space_type loc = *s;

    if( N > loc.capacity ) {

        if( loc.free_list != NULL ) {
            delete loc.free_list;
        }

        loc.free_list = new typename space_type::int_type [ free_list_size ];
        memset( loc.free_list, -1, free_list_size * sizeof( typename space_type::int_type ) );

        if( loc.free_map != NULL ) {
            delete loc.free_map;
        }

        loc.free_map = new unsigned int[ N ];
        memset( loc.free_map, 0, N * sizeof( unsigned int ) );
        
        loc.capacity = N;
    }
    loc.size = N;

    *s = loc;
}

template < class IntType, class OrderTag >
__global__ void _resize_space( device_free_space< IntType, OrderTag > * s, unsigned int N ) {
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    if( tid == 0 ) {
        _resize_space_impl( s, N );
    }
}

template < class IntType, class OrderTag >
__global__ void _delete_space( device_free_space< IntType, OrderTag > * s ) {
    if( s == NULL ) return;

    if( s->free_list != NULL ) {
        delete s->free_list;
        delete s->free_map;
    }
}

#endif  // DEVICE_FREE_SPACE_KERNELS_HPP_
