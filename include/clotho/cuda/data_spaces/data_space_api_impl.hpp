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
#ifndef DATA_SPACE_API_IMPL_HPP_
#define DATA_SPACE_API_IMPL_HPP_

#include "clotho/cuda/data_spaces/data_space_kernel_api.hpp"

template < class SpaceType >
void create_space( SpaceType *& space, unsigned int N ) {
    unsigned int n = sizeof( SpaceType );
    assert( cudaMalloc( (void **) &space, n ) == cudaSuccess );
    assert( cudaMemset( space, 0, n ) == cudaSuccess );

    if( N ) {
        resize_space( space, N );
    }
}

template < class SpaceType >
void resize_space( SpaceType * space, unsigned int N ) {
    _resize_space<<< 1, 1 >>>( space, N );
}

template < class SpaceType >
void delete_space( SpaceType * space ) {
    _delete_space<<< 1, 1 >>>( space );

    cudaFree( space );
}

#endif  // DATA_SPACE_API_IMPL_HPP_
