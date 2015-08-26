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
#ifndef POPCOUNT_KERNEL_H_
#define POPCOUNT_KERNEL_H_

#include <cuda.h>
#include "clotho/utility/popcount.hpp"

template < unsigned int W >
struct popcountGPU_constants;

template < >
struct popcountGPU_constants< 8 > : public popcount_constants< 8 > {
    typedef popcount_constants< 8 > base_type;
};

template < >
struct popcountGPU_constants< 4 > : public popcount_constants< 4 > {
    typedef popcount_constants< 4 > base_type;
};

template < class T >
struct popcountGPU : public popcountGPU_constants< sizeof( T ) > {
    typedef typename popcountGPU_constants< sizeof(T) >::base_type base_type;

    __device__ unsigned int evalGPU( T x ) {
        x -= (x >> 1) & base_type::M1;
        x = (x & base_type::M2) + ((x >> 2) & base_type::M2);
        x = (x + (x >> 4) ) & base_type::M4;
        return ( x * base_type::H01) >> base_type::DSHIFT;
    }
};

template < class T >
__global__ void computeHW( T * a, unsigned int N ) {
    int bid = blockIdx.y * gridDim.x + blockIdx.x;
    int tid = threadIdx.y * blockDim.x + threadIdx.x;

    int idx = bid * (blockDim.x * blockDim.y) + tid;

    popcountGPU< T > pc;

    if( idx < N ) {
        T x = a[idx];

        a[idx] = (T) pc.evalGPU(x);
    }
    __syncthreads();
}

#endif  // POPCOUNT_KERNEL_H_
