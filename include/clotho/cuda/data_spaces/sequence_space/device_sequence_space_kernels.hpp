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
#ifndef DEVICE_SEQUENCE_SPACE_KERNELS_HPP_
#define DEVICE_SEQUENCE_SPACE_KERNELS_HPP_

#include "clotho/cuda/data_spaces/sequence_space/device_sequence_space_def.hpp"
#include "clotho/cuda/data_spaces/sequence_space/device_sequence_space_kernel_api.hpp"

template < class IntType >
__device__ void _resize_space_impl( device_sequence_space< IntType > * sspace, unsigned int cols, unsigned int rows ) {
//    assert( blockIdx.y * gridDim.x + blockIdx.x == 0);
//    assert( threadIdx.y * blockDim.x + threadIdx.x == 0 );

    typedef typename device_sequence_space< IntType >::int_type int_type;
    unsigned int N = cols * rows;

    if( sspace->capacity < N ) {
        int_type * seqs = sspace->sequences;
        if( seqs ) {
            delete seqs;
        }

        seqs = new int_type[ N ];

        assert( seqs != NULL );

        sspace->sequences = seqs;
        sspace->capacity = N;
    }

    sspace->size = N;
    sspace->seq_count = rows;
    sspace->seq_width = cols;
}

template < class IntType >
__global__ void _resize_space( device_sequence_space< IntType > * sspace, unsigned int cols, unsigned int rows = 1 ) {
    assert( blockIdx.y * gridDim.x + blockIdx.x == 0 );

    if( threadIdx.y * blockDim.x + threadIdx.x == 0 ) {
        _resize_space_impl( sspace, cols, rows );
    }
};

template < class IntType, class ColumnSpaceType >
__global__ void _resize_space( device_sequence_space< IntType > * sspace, ColumnSpaceType * aspace, unsigned int seq_count ) {
    assert( blockIdx.y * gridDim.x + blockIdx.x == 0 );

    if( threadIdx.y * blockDim.x + threadIdx.x == 0 ) {
        typedef device_sequence_space< IntType > space_type;
        typedef typename space_type::int_type int_type;

        unsigned int W = aspace->capacity;
        W /= space_type::OBJECTS_PER_INT;

        _resize_space_impl( sspace, W, seq_count );
    }
}

template < class IntType >
__global__ void _delete_space( device_sequence_space< IntType > * sspace ) {
    if( sspace->sequences != NULL ) {
        delete sspace->sequences;
    }
}

#endif  // DEVICE_SEQUENCE_SPACE_KERNELS_HPP_
