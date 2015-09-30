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
#ifndef BASIC_DATA_SPACE_KERNELS_HPP_
#define BASIC_DATA_SPACE_KERNELS_HPP_

#include "clotho/cuda/data_spaces/basic_data_space.hpp"
#include "clotho/cuda/data_spaces/data_space_kernel_api.hpp"

template < class ValueType >
__global__ void _delete_space( basic_data_space< ValueType > * space ) {
    typename basic_data_space< ValueType >::value_type * tmp = space->data;

    if( tmp ) {
        delete tmp;
    }
}

template < class ValueType >
__device__ void _resize_space_impl( basic_data_space< ValueType > * space, unsigned int N ) {
    typedef basic_data_space< ValueType > space_type;

    if( space->capacity < N ) {
        if( space->data ) {
            delete space->data;
        }

        space->data = new typename space_type::value_type[ N ];
        memset( space->data, 0, N * sizeof( typename space_type::value_type) );
        space->capacity = N;
    }

    space->size = N;
}

template < class ValueType >
__global__ void _resize_space( basic_data_space< ValueType > * space, unsigned int N ) {
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    if( tid == 0 ) {
        _resize_space_impl( space, N );
    }
}

template < class ValueType, class SpaceType >
__global__ void _resize_space_to( basic_data_space< ValueType > * target, SpaceType * base ) {
    if( threadIdx.y * blockDim.x + threadIdx.x == 0 ) {
        unsigned int N = base->capacity;
        _resize_space_impl( target, N );
    }
}

#endif  // BASIC_DATA_SPACE_KERNELS_HPP_
