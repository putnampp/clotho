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
#ifndef ORDER_WARP_HPP_
#define ORDER_WARP_HPP_

#include <cuda.h>

#include "clotho/cuda/compute_capability.hpp"
#include "clotho/cuda/order_warp/order_warp_kernel.hpp"
#include "clotho/cuda/warp_sort.hpp"

namespace clotho {
namespace cuda {

/**
 * order_warp: produce periodically ordered key (& value) array(s).
 * Every WARP_SIZE elements of the key list (and their associated values)
 * are sorted into an ascending ordered and written to output vectors.
 *
 * This is effectively a wrapper object that invokes a CUDA GPU kernel
 * with appropriately sized execution configuration.
 */

template < class CC >
class order_warp;

template < >
struct order_warp< compute_capability< 3, 0 > > {
    typedef compute_capability< 3, 0 > comp_cap_type;

    order_warp() {}

/**
 *
 * All vectors are assumed to be allocated on the device.
 *
 * @param keys  1D vector of keys
 * @param values  1D vector of values
 * @param sorted_keys - 1D periodically ordered
 * @param sorted_values
 * @param N  size of the vectors
 *
 * keys and values sharing the same index are assumed to be a pair.
 */
    template < class K, class V >
    void sort( K * k, V * v, K * sorted_k, V * sorted_v, unsigned int N ) {
        while( N ) {
            // has full thread blocks
            dim3 blocks( 1, 1, 1), threads(1 ,1 , 1 );
            determine_configuration( N, blocks, threads );

            // invoke kernel
            unsigned int n = blocks.x * comp_cap_type::THREADS_PER_BLOCK;
            if( n > N ) {
                n = N;
            }

            comp_cap_type cc;
            order_warp_kernel<<< blocks, threads >>>( k, v, sorted_k, sorted_v, n, cc );

            k += n;
            v += n;
            sorted_k += n;
            sorted_v += n;
            N -= n;
        }
    }

    template < class K >
    void sort( K * k, K * sorted_k, unsigned int * sorted_v, unsigned int N ) {
        while( N ) {
            // has full thread blocks
            dim3 blocks( 1, 1, 1), threads(1 ,1 , 1 );
            determine_configuration( N, blocks, threads );

            // invoke kernel
            unsigned int n = blocks.x * comp_cap_type::THREADS_PER_BLOCK;
            if( n > N ) {
                n = N;
            }

            comp_cap_type cc;
            order_warp_kernel<<< blocks, threads >>>( k, sorted_k, sorted_v, n, cc );

            k += n;
            sorted_k += n;
            sorted_v += n;
            N -= n;
        }
    }
    void determine_configuration( size_t N, dim3 & blocks, dim3 & threads ) {
        blocks.x = N / comp_cap_type::THREADS_PER_BLOCK;
        if( blocks.x == 0 ) {
            blocks.x = 1;
        } else if( blocks.x > comp_cap_type::MAX_BLOCKS_X ) {
            blocks.x = comp_cap_type::MAX_BLOCKS_X;
        }
        
        threads.x = comp_cap_type::WARP_SIZE;
        threads.y = comp_cap_type::WARP_PER_BLOCK;
        threads.z = 1;
    }
};

}   // namespace cuda
}   // namespace clotho

#endif  // ORDER_WARP_HPP_
