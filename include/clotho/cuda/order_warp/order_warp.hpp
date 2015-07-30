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
 * order_warp: produce periodically ordered key & value arrays.
 * Every WARP_SIZE elements of the key list and their associated values
 * are ordered within the output vectors.
 *
 * Performs a bitonic sort using registers of threads within a warp.
 * The initial version of the algorithm is provided in:
 * 
 * @inproceedings{Demouth2013,
 * author = {Julien Demouth},
 * title = {Shuffle: Tips and Tricks},
 * year = {2013},
 * url = {http://on-demand.gputechconf.com/gtc/2013/presentations/S3174-Kepler-Shuffle-Tips-Tricks.pdf},
 * }
 *
 * @param keys  1D vector of keys [in]
 * @param values  1D vector of values
 * @param N  size of the vectors
 * @param sorted_keys - 1D periodically ordered
 * @param sorted_values
 *
 * keys and values sharing the same index are assumed to be a pair.
 */

template < class CC >
class order_warp;

template < >
struct order_warp< compute_capability< 3, 0 > > {
    typedef compute_capability< 3, 0 > comp_cap_type;

    order_warp() {}

    template < class K, class V >
    void sort( K * k, V * v, K * sorted_k, V * sorted_v, unsigned int N ) {
        unsigned int bcount = N / comp_cap_type::THREADS_PER_BLOCK;
        dim3 blocks( bcount, 1, 1), threads( comp_cap_type::WARP_SIZE, comp_cap_type::WARP_PER_BLOCK, 1 );

        comp_cap_type cc;
        if( bcount > 0 ) {
            // has full thread blocks

            // invoke kernel
            unsigned int n = bcount * comp_cap_type::THREADS_PER_BLOCK;
            order_warp_kernel<<< blocks, threads >>>( k, v, sorted_k, sorted_v, n, cc );

            k += n;
            v += n;
            sorted_k += n;
            sorted_v += n;
            N -= n;
        }

        if( N ) {
            // contains un-full thread block
            blocks.x = 1; blocks.y = 1; blocks.z = 1;
            order_warp_kernel<<< blocks, threads >>>( k, v, sorted_k, sorted_v, N, cc);
        }
    }
};

}   // namespace cuda
}   // namespace clotho

#endif  // ORDER_WARP_HPP_
