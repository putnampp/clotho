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
#ifndef WARP_KVSORT_HPP_
#define WARP_KVSORT_HPP_

#include <cuda.h>

#include "compute_capability.hpp"

namespace clotho {
namespace cuda {
/**
 * @param K - key
 * @param V - value
 * @param WW - warp width
 */
template < class K, class V, unsigned int WW >
struct warp_bitonic_sort;

// There is a known issue with sorting lists that are not a WARP_SIZE multiple
template < class K, class V >
struct warp_bitonic_sort< K, V, 32 > {
    typedef K key_type;
    typedef V value_type;
    typedef compute_capability< 3, 0 > comp_cap_type;

    static const unsigned int MAX_BITONIC_ROUNDS = 5;

    __device__ void sort( key_type & k, value_type & v, unsigned int round ) {
        unsigned int round_mask = ((round == MAX_BITONIC_ROUNDS) ? 0 : ((threadIdx.x >> round) & 0x01));
        unsigned int xor_mask = (1 << (round - 1));

        while( round-- ) {
            unsigned int dir = ((round_mask ^ (threadIdx.x >> round)) & 0x01);
            K sw_k = __shfl_xor( k, xor_mask );
            V sw_v = __shfl_xor( v, xor_mask );
            if( k < sw_k == dir ) {
                k = sw_k;
                v = sw_v;
            }

            xor_mask >>= 1;
        }
    }
};

template < class K, class V, class CC >
__global__ void kvsort( K * keys, V * values, K * sorted_keys, V * sorted_values, unsigned int N, CC comp_cap ) {

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int element_idx = (blockIdx.y * gridDim.x + blockIdx.x) * (blockDim.x * blockDim.y) + tid;

    typedef warp_bitonic_sort< K, V, CC::WARP_SIZE > sorter_type;
    sorter_type sorter;

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

    for( unsigned int i = 1; i <= sorter_type::MAX_BITONIC_ROUNDS; ++i ) {
    //    bitonic_sort( ke, va, i, element_idx, N );
        sorter.sort( ke, va, i);
    }

    if( element_idx < N ) {
        sorted_keys[ element_idx ] = ke;
        sorted_values[ element_idx ] = va;
    }
}

/**
 * warp_kvsort: produce periodically ordered key & value arrays.
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

template < class K, class V, class CC >
class warp_kvsort;

template < class K, class V >
struct warp_kvsort< K, V, compute_capability< 3, 0 > > {
    typedef K key_type;
    typedef V value_type;
    typedef compute_capability< 3, 0 > comp_cap_type;

//    static const unsigned int WARP_SIZE = compute_capability3::WARP_SIZE;
//    static const unsigned int THREADS_PER_BLOCK = ;
//    static const unsigned int WARPS_PER_BLOCK = 32;     // 1024 / 32
    static const unsigned int MAX_BITONIC_ROUNDS = 5;

    warp_kvsort() {}

    void sort( key_type * k, value_type * v, key_type * sorted_k, value_type * sorted_v, unsigned int N ) {
        unsigned int bcount = N / comp_cap_type::THREADS_PER_BLOCK;
        dim3 blocks( bcount, 1, 1), threads( comp_cap_type::WARP_SIZE, comp_cap_type::WARP_PER_BLOCK, 1 );

        comp_cap_type cc;
        if( bcount > 0 ) {
            // has full thread blocks

            // balance block dimensions

            // invoke kernel
            unsigned int n = bcount * comp_cap_type::THREADS_PER_BLOCK;
            kvsort<<< blocks, threads >>>( k, v, sorted_k, sorted_v, n, cc );

            k += n;
            v += n;
            sorted_k += n;
            sorted_v += n;
            N -= n;
        }

        if( N ) {
            // contains un-full thread block
            blocks.x = 1; blocks.y = 1; blocks.z = 1;
            kvsort<<< blocks, threads >>>( k, v, sorted_k, sorted_v, N, cc);
        }
    }
};

}   // namespace cuda
}   // namespace clotho

#endif  // WARP_KVSORT_HPP_
