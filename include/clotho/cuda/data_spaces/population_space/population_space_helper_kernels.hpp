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
#ifndef POPULATION_SPACE_HELPER_KERNELS_HPP_
#define POPULATION_SPACE_HELPER_KERNELS_HPP_

#include <cuda.h>

#include "clotho/cuda/data_spaces/sequence_space/device_sequence_space.hpp"
#include "clotho/cuda/data_spaces/free_space/device_free_space.hpp"

#include "clotho/cuda/data_spaces/tags/unordered_tag.hpp"
#include "clotho/cuda/data_spaces/tags/unit_ordered_tag.hpp"

#include "clotho/cuda/popcount_kernel.h"


template < class IntType >
__global__ void resize_free_space_kernel( device_sequence_space< IntType > * seq_space
                                        , device_free_space< IntType, unordered_tag > * free_space ) {
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    unsigned int seq_width = seq_space->seq_width;
    if( tid == 0 && seq_width != 0 ) {
        _resize_space_impl( free_space, seq_width * device_sequence_space< IntType >::OBJECTS_PER_INT );
    }
}
/*
template < class IntType >
__global__ void update_free_space_kernel( device_sequence_space< IntType > * seq_space
                                        , device_free_space< IntType, unordered_tag > * free_space ) {
    typedef device_sequence_space< IntType >        sequence_space_type;
    typedef typename sequence_space_type::int_type  sequence_int_type;

    typedef device_free_space< IntType, unordered_tag >  free_space_type;
    typedef typename free_space_type::int_type      free_int_type;

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    unsigned int seq_width = seq_space->seq_width;
    unsigned int seq_count = seq_space->seq_count;

    if( seq_count == 0 ) return;

    unsigned int column = tid;

    free_int_type * flist = free_space->free_list;

    while( column < seq_width ) {
        sequence_int_type * seq = seq_space->sequences;

        sequence_int_type union_mask = seq[column];
        sequence_int_type int_mask = union_mask;

        unsigned int i = 1;
        while( i < seq_count ) {
            sequence_int_type s = *seq;

            union_mask |= s;
            int_mask &= s;

            seq += seq_width;
            ++i;
        }

        flist[column] = (int_mask | (~union_mask));

        column += blockDim.x;
    }
}*/

template < class IntType >
__global__ void update_free_space_kernel2( device_sequence_space< IntType > * seq_space
                                        , device_free_space< IntType, unordered_tag > * free_space ) {

    typedef device_sequence_space< IntType >        sequence_space_type;
    typedef typename sequence_space_type::int_type  sequence_int_type;

    typedef device_free_space< IntType, unordered_tag >  free_space_type;
    typedef typename free_space_type::int_type      free_int_type;

    unsigned int seq_width = seq_space->seq_width;
    unsigned int seq_count = seq_space->seq_count;

    if( seq_count == 0 || seq_width == 0 ) return;

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;

    unsigned int column = bid;

    free_int_type * flist = free_space->free_list;
    free_int_type * xlist = free_space->fixed_list;
    free_int_type * llist = free_space->lost_list;

    sequence_int_type * seqs = seq_space->sequences;

    unsigned int row_offset = (blockDim.x * blockDim.y) * seq_width;
    unsigned int col_offset = (gridDim.x * gridDim.y);

    unsigned int _count = seq_space->size;

    while( column < seq_width ) {
        unsigned int ridx = tid * seq_width + column;

        free_int_type union_mask = 0, inter_mask = (free_int_type) -1;
        if( ridx < _count ) {
            sequence_int_type v = seqs[ ridx ];
            union_mask = v;
            inter_mask = v;
            ridx += row_offset;

            while( ridx < _count ) {
                v = seqs[ ridx ];
                union_mask |= v;
                inter_mask &= v;
                ridx += row_offset;
            }
        }
        __syncthreads();

        // reduce the masks
        for( unsigned int i = 1; i < 32; i <<= 1 ) {
            sequence_int_type t = __shfl_up( inter_mask, i);
            sequence_int_type u = __shfl_up( union_mask, i);

            inter_mask &= ((tid >= i) ? t : (sequence_int_type) -1);
            union_mask |= ((tid >= i) ? u : 0 );
        }

        if( tid == 31 ) {
            llist[column] = ~union_mask;
            xlist[column] = inter_mask;
            flist[column] = (inter_mask | (~union_mask));
        }
        __syncthreads();
        column += col_offset;
    }
}

/**
 *
 * Intended to be executed two blocks, single warp
 */
template < class IntType, class OrderTag >
__global__ void update_free_space_total_kernel( device_free_space< IntType, OrderTag > * fspace ) {
    typedef device_free_space< IntType, OrderTag >      space_type;
    typedef typename space_type::int_type               int_type;

    unsigned int bid = blockIdx.y * gridDim.x + blockDim.x;

    if( bid >= 2 ) return;
    popcountGPU< int_type > pc;

    unsigned int count = 0;
    unsigned int width = fspace->size;

    width /= (sizeof( int_type ) * 8);

    unsigned int idx = threadIdx.y * blockDim.x + threadIdx.x;

    int_type * _list;
    if(bid == 0 ) { // true for all threads in block
        _list = fspace->free_list;
    } else {
        _list = fspace->fixed_list;
    }
    __syncthreads();

    while( idx < width ) {
        int_type c = _list[idx];

        count += pc.evalGPU( c );

        idx += (blockDim.x * blockDim.y);
    }
    __syncthreads();

    idx = threadIdx.y * blockDim.x + threadIdx.x;

    for( unsigned int i = 1; i < 32; i <<= 1 ) {
        unsigned int c = __shfl_up( count, i );
        count += ( idx >= i ) * c;
    }

    if( idx == 31 ) {
        if( bid == 0 ) {
            fspace->total = count;
        } else {
            fspace->fixed_count = count;
        }
    }
}

template < class IntType >
__global__ void update_free_map_kernel( device_free_space< IntType, unordered_tag > * free_space ) {
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    unsigned int warp_id = (tid >> 5);
    unsigned int lane_id = (tid & 31);

    typedef device_free_space< IntType, unordered_tag > space_type;
    typedef typename space_type::int_type               int_type;

    __shared__ unsigned int s_map[ 32 ];

    space_type loc = *free_space;

    popcountGPU< int_type > pc;

    unsigned int count = 0, prev_sum = 0;
    unsigned int i = warp_id;

    unsigned int N = (loc.size / space_type::OBJECTS_PER_INT);

    unsigned int padded_N = ((N & 31) ? ((N & (~31)) + 32) : N);

    while( i < padded_N ) {
        int_type _c = 0;
        if( i < N ) {
            _c = loc.free_list[i];
        }
        __syncthreads();

        count += pc.evalGPU( _c );

        unsigned int offset = ((_c >> lane_id) & 1);

        bool is_free = (offset != 0);

        for( unsigned int j = 1; j < 32; j <<= 1 ) {
            unsigned int t = __shfl_up( offset, j );
            offset += ( lane_id >= j ) * t;
        }

        if( lane_id == 31 ) {
            s_map[ warp_id ] = offset;
        }
        __syncthreads();

        unsigned int tmp = s_map[lane_id];

        for( unsigned int j = 1; j < 32; j <<= 1 ) {
            unsigned int t = __shfl_up( tmp, j );
            tmp += ( lane_id >= j ) * t;
        }

        if( warp_id > 0 ) {
            unsigned int t = __shfl( tmp, warp_id - 1);
            offset += t;
        }
        __syncthreads();

        offset += prev_sum;

        if( is_free ) {
            loc.free_map[offset - 1] = i * 32 + lane_id;
        }
        __syncthreads();

        unsigned int t = __shfl( tmp, 31 );
        prev_sum = t;
        i += 32;
    }

    if( lane_id == 0 ) {
        s_map[ warp_id ] = count;
    }
    __syncthreads();

    if( warp_id == 0 ) {
        count = s_map[ lane_id ];

        for( i = 1; i < 32; i <<= 1 ) {
            unsigned int t = __shfl_up( count, i );
            count += (tid >= i) * t;
        }

        if( tid == 31 ) {
            free_space->total = count;
        }
    }
}

/*
template < class IntType, class OrderTag >
__global__ void clear_free_space_kernel( device_sequence_space< IntType > * seq_space, device_free_space< IntType, OrderTag > * free_space ) {
    typedef device_sequence_space< IntType >        sequence_space_type;
    typedef typename sequence_space_type::int_type  sequence_int_type;

    typedef device_free_space< IntType, OrderTag >  free_space_type;
    typedef typename free_space_type::int_type      free_int_type;

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    unsigned int seq_width = seq_space->seq_width;
    unsigned int seq_count = seq_space->seq_count;

    if( seq_count == 0 ) return;

    unsigned int column = tid;

    free_int_type * flist = free_space->free_list + tid;

    while( column < seq_width ) {

        sequence_int_type * seq = seq_space->sequences + column;

        free_int_type mask = *flist;

        mask = ~mask;   // keep all variable (ie. not fixed or lost)

        unsigned int i = 0;
        while( i < seq_count ) {
            sequence_int_type s = *seq;

            s &= mask;

            *seq = s;

            seq += seq_width;
            ++i;
        }

        flist += blockDim.x;
        column += blockDim.x;
    }
}*/

/*
template < class IntType, class OrderTag >
__global__ void clear_free_space_kernel2( device_sequence_space< IntType > * seq_space, device_free_space< IntType, OrderTag > * free_space ) {
    typedef device_sequence_space< IntType >        sequence_space_type;
    typedef typename sequence_space_type::int_type  sequence_int_type;

    typedef device_free_space< IntType, OrderTag >  free_space_type;
    typedef typename free_space_type::int_type      free_int_type;

    unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;

    unsigned int seq_count = seq_space->seq_count;
    if( bid >= seq_count ) { return; }

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int seq_width = seq_space->seq_width;

    unsigned int seq_offset = gridDim.x * gridDim.y;
    unsigned int flist_offset = blockDim.x * blockDim.y;

    sequence_int_type * seqs = seq_space->sequences;
    free_int_type * flist = free_space->free_list;

    while( bid < seq_count ) {

        unsigned int sidx = bid * seq_width + tid;
        unsigned int j = tid;
        while( j < seq_width ) {
            free_int_type mask = flist[j];
            sequence_int_type s = seqs[sidx];

            seqs[sidx] = (s & ~mask);

            j += flist_offset;
            sidx += flist_offset;
        }
        __syncthreads();
        bid += seq_offset;
    }
}*/

template < class IntType, class OrderTag >
__global__ void clear_free_space_kernel3( device_sequence_space< IntType > * seq_space, device_free_space< IntType, OrderTag > * free_space ) {
    typedef device_sequence_space< IntType >        sequence_space_type;
    typedef typename sequence_space_type::int_type  sequence_int_type;

    typedef device_free_space< IntType, OrderTag >  free_space_type;
    typedef typename free_space_type::int_type      free_int_type;

    unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;
    unsigned int seq_width = seq_space->seq_width;

    if( bid >= seq_width ) { return; }  // get out quick for unnecessary blocks

    unsigned int _size = seq_space->size;
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    unsigned int seq_offset = blockDim.x * blockDim.y * seq_width;
    unsigned int flist_offset = gridDim.x * gridDim.y;

    sequence_int_type * seqs = seq_space->sequences;
    free_int_type * flist = free_space->free_list;

    while( bid < seq_width ) {  // true for all threads in block
        free_int_type mask = flist[ bid ];
        mask = ~mask;

        unsigned int sidx = tid * seq_width + bid;
        while( sidx < _size ) { // threads access memory in different strides
            sequence_int_type s = seqs[sidx];

            s &= mask;
            seqs[sidx] = s;

            sidx += seq_offset;
        }
        __syncthreads();

        bid += flist_offset;
    }
}

#endif  // POPULATION_SPACE_HELPER_KERNELS_HPP_