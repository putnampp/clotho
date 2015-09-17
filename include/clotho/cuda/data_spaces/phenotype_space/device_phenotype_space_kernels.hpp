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
#ifndef DEVICE_PHENOTYPE_SPACE_KERNELS_HPP_
#define DEVICE_PHENOTYPE_SPACE_KERNELS_HPP_

#include "clotho/cuda/data_spaces/phenotype_space/device_phenotype_space_def.hpp"
#include "clotho/cuda/data_spaces/data_space_kernel_api.hpp"

template < class RealType >
__global__ void _delete_space( device_phenotype_space< RealType > * space ) {
    typename device_phenotype_space< RealType >::real_type * tmp = space->data;

    if( tmp ) {
        delete tmp;
    }
}

/*template < class RealType >
__device__ void _resize_space_impl( device_phenotype_space< RealType > * space, unsigned int N ) {
    typedef device_phenotype_space< RealType > space_type;

    if( space->capacity < N ) {
        typename space_type::real_type * tmp = space->data;
        if( tmp ) {
            delete tmp;
        }

        space->data = new typename space_type::real_type[ N ];
        memset( space->data, 0, N * sizeof( typename space_type::real_type) );
        space->capacity = N;
    }

    space->size = N;
}*/

template < class RealType >
__global__ void _resize_space( device_phenotype_space< RealType > * space, unsigned int N ) {
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    if( tid == 0 ) {
        _resize_space_impl( (basic_data_space< RealType > *) space, N );
    }
}

#endif  // DEVICE_PHENOTYPE_SPACE_KERNELS_HPP_
