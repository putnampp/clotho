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
#ifndef QUADRATIC_FITNESS_KERNEL_HPP_
#define QUADRATIC_FITNESS_KERNEL_HPP_

#include "clotho/cuda/data_spaces/basic_data_space.hpp"
#include <cuda.h>

template < class RealType >
__global__ void quadratic_fitness_kernel( basic_data_space< RealType > * phenos, RealType scale_coeff, basic_data_space< RealType > * fitness ) {
    typedef RealType real_type;

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    unsigned int N = phenos->size;
    unsigned int M = fitness->size;

    assert( N <= M );

    real_type * pdata = phenos->data;
    real_type * fdata = fitness->data;

    while( tid < N ) {
        real_type x = pdata[tid];

        x /= scale_coeff;
        x *= x;

        // x = 1.0 - ((real_type)(x > 1.0)* 1.0) - ((real_type)(x < 1.0) * x)
        x = (( x > 1.0 ) ? 0.0 : (1.0 - x));

        fdata[tid] = x;
        tid += (blockDim.x * blockDim.y);
    }
}

#endif  // QUADRATIC_FITNESS_KERNEL_HPP_
