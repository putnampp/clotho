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
#ifndef ORDER_WARP_KERNEL_HPP_
#define ORDER_WARP_KERNEL_HPP_

#include "clotho/cuda/warp_sort.hpp"

namespace clotho {
namespace cuda {

template < class K, class V, class CC >
__global__ void order_warp_kernel( K * keys, V * values, K * sorted_keys, V * sorted_values, unsigned int N, CC comp_cap ) {
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int element_idx = (blockIdx.y * gridDim.x + blockIdx.x) * (blockDim.x * blockDim.y) + tid;

    K ke = 0;
    V va = 0;
    if( element_idx < N ) {
        ke = keys[ element_idx ];
        va = values[ element_idx ];
    }
    __syncthreads();

    if( N % CC::WARP_SIZE ) {
        K max = ke;
        V vmax = va;

        for( unsigned int i = 1; i < CC::WARP_SIZE; i <<= 1 ) {
            K tmp = __shfl_down( max, i );
            V tmpv = __shfl_down( vmax, i );
            if( ((element_idx + i) < N) && (tmp > max) ) {
                max = tmp;
                vmax = tmpv;
            }
        }

        K tmp = __shfl( max, 0);
        V tmpv = __shfl( vmax, 0);
        if( element_idx >= N ) {
            ke = tmp;
            va = tmpv;
        }
    }
    __syncthreads();

    typedef warp_sort< CC::WARP_SIZE > sorter_type;
    sorter_type sorter;
    for( unsigned int i = 1; i <= sorter_type::MAX_BITONIC_ROUNDS; ++i ) {
        sorter.sort( ke, va, i);
    }

    if( element_idx < N ) {
        sorted_keys[ element_idx ] = ke;
        sorted_values[ element_idx ] = va;
    }
}

}   // namespace cuda
}   // namespace clotho

#endif  // ORDER_WARP_KERNEL_HPP_
