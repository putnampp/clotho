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
#ifndef SEQUENCE_HAMMING_WEIGHT_KERNEL_HPP_
#define SEQUENCE_HAMMING_WEIGHT_KERNEL_HPP_

#include <cuda.h>

#include "clotho/cuda/data_spaces/sequence_space/device_sequence_space.hpp"
#include "clotho/cuda/data_spaces/basic_data_space.hpp"

#include "clotho/cuda/popcount_kernel.h"
#include "clotho/utility/algorithm_version.hpp"

template < class IntType >
__global__ void _resize_space_to( basic_data_space< unsigned int > * res, device_sequence_space< IntType > * seqs ) {
    if( threadIdx.y * blockDim.x + threadIdx.x == 0 ) {
        unsigned int N = seqs->seq_count;
        _resize_space_impl( res, N );
    }
}

template < class IntType >
__global__ void sequence_hamming_weight_kernel( device_sequence_space< IntType > * seqs, basic_data_space< unsigned int > * res, clotho::utility::algo_version< 1 > * v ) {

    typedef device_sequence_space< IntType >    space_type;
    typedef typename space_type::int_type       int_type;

    unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;
    unsigned int _count = seqs->seq_count;

    if( bid >= _count ) return;

    assert( _count <= res->size );  // sanity check: enough allocated space for results

    popcountGPU< int_type > pc;

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int tpb = blockDim.x * blockDim.y; // threads per block
    unsigned int bpg = gridDim.x * gridDim.y;   // blocks per grid

    assert( (tpb & 31) == 0); // sanity check: all warps are full

    unsigned int lane_id = (tid & 31);
    unsigned int warp_id = (tid >> 5);

    unsigned int _width = seqs->seq_width;

    int_type * sptr = seqs->sequences;

    unsigned int * countptr = res->data;

    __shared__ unsigned int buffer[32];

    while( bid < _count ) {
        unsigned int degree = 0;

        unsigned int end = (bid + 1) * _width;
        unsigned int idx = (end - _width) + tid;

        while( idx < end ) {    // at most 1 warp diverges
            int_type b = sptr[ idx ];
            degree += pc.evalGPU( b );

            idx += tpb;
        }
        __syncthreads();

        for( unsigned int i = 1; i < 32; i <<= 1 ) {
            unsigned int d = __shfl_up(degree, i );
            degree += (( lane_id >= i ) * d);
        }

        if( tpb > 32 ) { // true for all threads in block
            if( lane_id == 31 ) {
                buffer[ warp_id ] = degree;
            }
            __syncthreads();

            degree = buffer[ lane_id ];
            __syncthreads();

            for( unsigned int i = 1; i < 32; i <<= 1 ) {
                unsigned int d = __shfl_up( degree, i );
                degree +=  ((lane_id >= i) * d);
            }
        }

        if( tid == 31 ) {
            countptr[ bid ] = degree;
        }
        __syncthreads();
        bid += bpg;
    }
}

template < class IntType >
__global__ void sequence_hamming_weight_kernel( device_sequence_space< IntType > * seqs, basic_data_space< unsigned int > * res, clotho::utility::algo_version< 2 > * v ) {

    typedef device_sequence_space< IntType >    space_type;
    typedef typename space_type::int_type       int_type;

    unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;
    unsigned int _count = seqs->seq_count;

    if( bid >= _count ) return;

    assert( _count <= res->size );  // sanity check: enough allocated space for results

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int tpb = blockDim.x * blockDim.y; // threads per block
    unsigned int bpg = gridDim.x * gridDim.y;   // blocks per grid
    unsigned int wpb = (tpb >> 5);

    assert( (tpb & 31) == 0); // sanity check: all warps are full

    unsigned int lane_id = (tid & 31);
    unsigned int warp_id = (tid >> 5);

    unsigned int _width = seqs->seq_width;

    int_type * sptr = seqs->sequences;

    unsigned int * countptr = res->data;

    __shared__ unsigned int buffer[32];

    while( bid < _count ) {
        unsigned int degree = 0;

        unsigned int end = (bid + 1) * _width;
        unsigned int idx = (end - _width) + warp_id;

        while( idx < end ) {    // all threads in a warp read same int; no divergence
            int_type b = sptr[ idx ];
            degree += ((b >> lane_id) & 1);
            idx += wpb;
        }
        __syncthreads();

        for( unsigned int i = 1; i < 32; i <<= 1 ) {
            unsigned int d = __shfl_up(degree, i );
            degree += (( lane_id >= i ) * d);
        }

        if( tpb > 32 ) { // true for all threads in block
            if( lane_id == 31 ) {
                buffer[ warp_id ] = degree;
            }
            __syncthreads();

            degree = buffer[ lane_id ];
            __syncthreads();

            for( unsigned int i = 1; i < 32; i <<= 1 ) {
                unsigned int d = __shfl_up( degree, i );
                degree +=  ((lane_id >= i) * d);
            }
        }

        if( tid == 31 ) {
            countptr[ bid ] = degree;
        }
        __syncthreads();
        bid += bpg;
    }
}

template < class IntType >
__global__ void sequence_hamming_weight_kernel( device_sequence_space< IntType > * seqs, basic_data_space< unsigned int > * res, clotho::utility::algo_version< 3 > * v ) {

    typedef device_sequence_space< IntType >    space_type;
    typedef typename space_type::int_type       int_type;

    unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;
    unsigned int _count = seqs->seq_count;

    if( bid >= _count ) return;

    assert( _count <= res->size );  // sanity check: enough allocated space for results

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int tpb = blockDim.x * blockDim.y; // threads per block
    unsigned int bpg = gridDim.x * gridDim.y;   // blocks per grid
    unsigned int wpb = (tpb >> 5);
    unsigned int spg = wpb * bpg;   // sequences/grid = sequences(warps)/block * blocks/grid

    assert( (tpb & 31) == 0); // sanity check: all warps are full

    unsigned int lane_id = (tid & 31);
    unsigned int warp_id = (tid >> 5);

    unsigned int _width = seqs->seq_width;

    int_type * sptr = seqs->sequences;

    unsigned int * countptr = res->data;

    unsigned int max_seq_id = _count / wpb; // max_rounds = sequences * block/sequences 
    max_seq_id += ((_count % wpb) ? 1 : 0); // would !!(_count % spg) be more efficient?
    max_seq_id *= wpb;

    unsigned int seq_id = bid * wpb + warp_id;

    while( seq_id < max_seq_id ) {  // blocks of grid may terminate early; only block for tail may diverge
        unsigned int degree = 0;

        unsigned int end = (seq_id + 1) * _width;
        unsigned int idx = end - ((seq_id < _count) ? _width : 0);  // true for all threads in warp

        while( idx < end ) {    // all threads in a warp read same bit block; no divergence
            int_type b = sptr[ idx++ ]; // all threads in a warp read/load same bit block
            degree += ((b >> lane_id) & 1); // would !!( b & lane_mask), where (lane_mask = (1 << lane_id)), be more efficient?
        }
        __syncthreads();    // sync all warps

        for( unsigned int i = 1; i < 32; i <<= 1 ) {
            unsigned int d = __shfl_up(degree, i );
            degree += (( lane_id >= i ) * d);
        }

        if( lane_id == 31 && (seq_id < _count) ) {
            countptr[ seq_id ] = degree;
        }
        __syncthreads();

        seq_id += spg;
    }
}

template < class IntType >
__global__ void sequence_hamming_weight_kernel( device_sequence_space< IntType > * seqs, basic_data_space< unsigned int > * res, clotho::utility::algo_version< 4 > * v ) {

    typedef device_sequence_space< IntType >    space_type;
    typedef typename space_type::int_type       int_type;

    unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;
    unsigned int _count = seqs->seq_count;

    if( bid >= _count ) return;

    assert( _count <= res->size );  // sanity check: enough allocated space for results

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int tpb = blockDim.x * blockDim.y; // threads per block
    unsigned int bpg = gridDim.x * gridDim.y;   // blocks per grid
    unsigned int wpb = (tpb >> 5);
    unsigned int spg = wpb * bpg;   // sequences/grid = sequences(warps)/block * blocks/grid

    assert( (tpb & 31) == 0); // sanity check: all warps are full

    unsigned int lane_id = (tid & 31);
    unsigned int warp_id = (tid >> 5);

    unsigned int _width = seqs->seq_width;

    int_type * sptr = seqs->sequences;

    unsigned int * countptr = res->data;

    unsigned int max_seq_id = _count / wpb; // max_rounds = sequences * block/sequences 
    max_seq_id += ((_count % wpb) ? 1 : 0); // would !!(_count % wpb) be more efficient?
    max_seq_id *= wpb;

    unsigned int seq_id = bid * wpb + warp_id;

    popcountGPU< int_type > pc;

    while( seq_id < max_seq_id ) {  // blocks of grid may terminate early; only block for tail may diverge
        unsigned int degree = 0;

        unsigned int end = (seq_id + 1) * _width;
        unsigned int idx = end - ((seq_id < _count) ? _width : 0);  // true for all threads in warp
        idx += lane_id;

        while( idx < end ) {    // allows thread divergence
            int_type b = sptr[ idx ]; // all threads in a warp read/load same bit block
            degree += pc.evalGPU( b );
            idx += 32;  // warp_size
        }
        __syncthreads();    // sync all warps

        for( unsigned int i = 1; i < 32; i <<= 1 ) {
            unsigned int d = __shfl_up(degree, i );
            degree += (( lane_id >= i ) * d);
        }

        if( lane_id == 31 && (seq_id < _count) ) {
            countptr[ seq_id ] = degree;
        }
        __syncthreads();

        seq_id += spg;
    }
}

/**
 *
 * 1 block per sequence
 *
 */
template < class IntType >
__global__ void sequence_hamming_weight_kernel( device_sequence_space< IntType > * seqs, basic_data_space< unsigned int > * res, clotho::utility::algo_version< 5 > * v ) {
    assert( blockDim.x == 32 );

    __shared__ unsigned int buffer[ 32 ];

    if( threadIdx.y == 0 ) {
        buffer[ threadIdx.x ] = 0;
    }
    __syncthreads();

    const unsigned int WIDTH = seqs->seq_width;

    IntType * seq_ptr = seqs->sequences;

    unsigned int N = 0;
    unsigned int seq_idx = blockIdx.y * gridDim.x  + blockIdx.x;
    unsigned int seq_begin = seq_idx * WIDTH + threadIdx.y;
    unsigned int seq_end = (seq_idx + 1) * WIDTH;
    while( seq_begin < seq_end ) {
        IntType x = seq_ptr[ seq_begin ];
        N += (( x >> threadIdx.x) & 1);
        seq_begin += blockDim.y;    
    }
    __syncthreads();

    for( unsigned int i = 1; i < 32; i <<= 1 ) {
        unsigned int t = __shfl_up( N, i );
        N += ((unsigned int) (threadIdx.x >= i) * t);
    }

    
    if( threadIdx.x == 31 ) {
        buffer[ threadIdx.y ] = N;
    }
    __syncthreads();

    if( threadIdx.y == 0 ) {
        N = buffer[ threadIdx.x ];

        for( unsigned int i = 0; i < 32; i <<= 1 ) {
            unsigned int t = __shfl_up( N, i );
            N += ((unsigned int) (threadIdx.x >= i) * t);
        }

        if( threadIdx.x == 31 ) {
            res->data[ seq_idx ] = N;           
        }
    }
}

#endif  // SEQUENCE_HAMMING_WEIGHT_KERNEL_HPP_
