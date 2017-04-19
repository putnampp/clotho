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
#ifndef DEVICE_FREE_SPACE_KERNELS_HPP_
#define DEVICE_FREE_SPACE_KERNELS_HPP_

#include "clotho/cuda/data_spaces/free_space/device_free_space_def.hpp"
#include "clotho/cuda/data_spaces/data_space_kernel_api.hpp"

template < class IntType, class OrderTag >
__device__ void _resize_space_impl( device_free_space< IntType, OrderTag > * s, unsigned int N ) {
//    printf("Resize free space: %d\n", N );
    assert( blockIdx.y * gridDim.x + blockIdx.x == 0 );
    assert( threadIdx.y * blockDim.x + threadIdx.x == 0 );

    typedef device_free_space< IntType, OrderTag > space_type;
    typedef typename space_type::int_type   int_type;

    unsigned int free_list_size = N / space_type::OBJECTS_PER_INT;
    if( N % space_type::OBJECTS_PER_INT ) {
        ++free_list_size;
        N = free_list_size * space_type::OBJECTS_PER_INT;
    }

    if( N > s->capacity ) {

        if( s->free_list != NULL ) {
            delete [] s->free_list;
            delete [] s->fixed_list;
            delete [] s->lost_list;
            delete [] s->free_map;
        }

        int_type * flist = new int_type [ free_list_size ];
        memset( flist, -1, free_list_size * sizeof( int_type ) );

        int_type * xlist = new int_type[ free_list_size ];
        memset( xlist, 0, free_list_size * sizeof( int_type ) );

        int_type * llist = new int_type[ free_list_size ];
        memset( llist, 0, free_list_size * sizeof( int_type ) );

        unsigned int * fmap = new unsigned int[ N ];
        memset( fmap, 0, N * sizeof( unsigned int ) );
        
        s->free_list = flist;
        s->fixed_list = xlist;
        s->lost_list = llist;

        s->free_map = fmap;
        s->capacity = N;
    }
    s->size = N;
}

template < class IntType, class OrderTag >
__global__ void _resize_space( device_free_space< IntType, OrderTag > * s, unsigned int N ) {
    assert( blockIdx.y * gridDim.x + blockIdx.x == 0 );

    if( threadIdx.y * blockDim.x + threadIdx.x == 0 ) {
        _resize_space_impl( s, N );
    }
}

template < class IntType, class OrderTag >
__global__ void _delete_space( device_free_space< IntType, OrderTag > * s ) {
    assert ( blockIdx.y * gridDim.x + blockIdx.x == 0 );
    assert ( threadIdx.y * blockDim.x + threadIdx.x == 0 );

    if( s == NULL ) return;

    if( s->free_list != NULL ) {
        delete [] s->free_list;
        delete [] s->fixed_list;
        delete [] s->lost_list;
        delete [] s->free_map;
    }

    
    s->free_list = NULL;
    s->fixed_list = NULL;
    s->lost_list = NULL;
    s->free_map = NULL;
}

template < class IntType, class OrderTag >
__device__ void _update_space( device_free_space< IntType, OrderTag > * in_space
                                , device_free_space< IntType, OrderTag > * out_space ) {

    typedef device_free_space< IntType, OrderTag > space_type;
    typedef typename space_type::int_type int_type;

    unsigned int N = in_space->capacity;
    unsigned int M = out_space->capacity;

    assert( N <= M );

    N /= (sizeof(int_type) * 8);
    N *= sizeof(int_type);

    M /= (sizeof(int_type) * 8);
    M *= sizeof(int_type);

    if( N == 0 || M == 0 ) {
        //printf( "Unexpected update of free space: M= %d; N=%d\n", M, N );
        return;
    }

    //printf("updating free space: %d -> %d\n", N, M );
    // reset out_space free list
    memset( out_space->free_list, 255, M );
    memset( out_space->fixed_list, 0, M );
    memset( out_space->lost_list, 0, M );

    // copy the in_space free list into out_space
    memcpy( out_space->free_list, in_space->free_list, N );
    memcpy( out_space->fixed_list, in_space->fixed_list, N );
    memcpy( out_space->lost_list, in_space->lost_list, N );
}

template < class IntType, class OrderTag >
__global__ void update_space_kernel( device_free_space< IntType, OrderTag > * in_space
                                        , device_free_space< IntType, OrderTag > * out_space ) {

    typedef device_free_space< IntType, OrderTag > space_type;
    typedef typename space_type::int_type int_type;

    unsigned int N = in_space->capacity;
    unsigned int M = out_space->capacity;

    assert( N <= M );

    assert( N % (sizeof(int_type) * 8) == 0 );
    assert( M % (sizeof(int_type) * 8) == 0 );

    N /= (sizeof(int_type) * 8);
    M /= (sizeof(int_type) * 8);

    if( N == 0 || M == 0 ) {
        //printf( "Unexpected update of free space: M= %d; N=%d\n", M, N );
        return;
    }

    int_type * in_data, *out_data;
    unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;
    if( bid == 0 ) {    // true for all threads
        in_data = in_space->free_list;
        out_data = out_space->free_list;
    } else if ( bid == 1 ) {    // true for all threads
        in_data = in_space->fixed_list;
        out_data = out_space->fixed_list;
    } else if ( bid == 2 ) {    // true for all threads
        in_data = in_space->lost_list;
        out_data = out_space->lost_list;
    } else {    // true for all threads
        return;
    }

    unsigned int tpb = blockDim.x * blockDim.y;
    unsigned int idx = threadIdx.y * blockDim.x + threadIdx.x;
    while( idx < N ) {
        out_data[idx] = in_data[idx];
        idx += tpb;
    }

    int_type v = 0;
    if( bid == 0 ) { v = (~v); }
    __syncthreads();

    while( idx < M ) {
        out_data[idx] = v;
        idx += tpb;
    }
}
#endif  // DEVICE_FREE_SPACE_KERNELS_HPP_
