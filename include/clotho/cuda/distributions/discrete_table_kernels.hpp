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
#ifndef DISCRETE_TABLE_KERNELS_HPP_
#define DISCRETE_TABLE_KERNELS_HPP_

#include "clotho/cuda/distributions/discrete_table.hpp"
#include "clotho/cuda/data_spaces/data_space_kernel_api.hpp"

template < class IntType, class RealType >
__device__ void _resize_space_impl( discrete_table< IntType, RealType > * tbl, unsigned int N ) {
    typedef discrete_table< IntType, RealType > table_type;
    typedef typename table_type::real_type      real_type;
    typedef typename table_type::int_type       int_type;
    
    if( N > tbl->capacity ) {
        if( tbl->threshold != NULL ) {
            delete tbl->threshold;
            delete tbl->alternative;
        }

        tbl->threshold = new real_type[ N ];
        assert( tbl->threshold != NULL );

        tbl->alternative = new int_type[ N ];
        assert( tbl->alternative != NULL );

        tbl->capacity = N;
    }
    tbl->size = N;
}

template < class IntType, class RealType >
__global__ void _resize_space( discrete_table< IntType, RealType > * tbl, unsigned int N ) {
    assert( blockIdx.y * gridDim.x + blockIdx.x == 0 );

    if( threadIdx.y * blockDim.x + threadIdx.x == 0 ) {
        _resize_space_impl( tbl, N );
    }
}

template < class IntType, class RealType >
__device__ void _delete_space_impl( discrete_table< IntType, RealType > * tbl ) {
    if( blockIdx.y * gridDim.x + blockIdx.x != 0 ) return;
    if( threadIdx.y * blockDim.x + threadIdx.x != 0 ) return;

    if( tbl->threshold != NULL ) {
        delete tbl->threshold;
        delete tbl->alternative;
    }

    tbl->threshold = NULL;
    tbl->alternative = NULL;
    tbl->capacity = 0;
    tbl->size = 0;
}

#endif  // DISCRETE_TABLE_KERNELS_HPP_
