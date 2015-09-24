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
__global__ void _resize_space( device_sequence_space< IntType > * sspace, unsigned int N ) {
    if( sspace->capacity < N ) {
        if( sspace->sequences ) {
            delete sspace->sequences;
        }

        sspace->sequences = new typename device_sequence_space< IntType >::int_type[ N ];

        memset( sspace->sequences, 0, N * sizeof( typename device_sequence_space< IntType >::int_type ) );
        sspace->capacity = N;
    }
    sspace->size = N;
};

template < class IntType, class ColumnSpaceType >
__global__ void _resize_space( device_sequence_space< IntType > * sspace, ColumnSpaceType * aspace, unsigned int seq_count ) {
    typedef device_sequence_space< IntType > space_type;
    typedef typename space_type::int_type int_type;

    unsigned int W = aspace->capacity;
    W /= space_type::OBJECTS_PER_INT;
    unsigned int N = seq_count * W;

    if( sspace->capacity < N ) {
        //printf("Resizing sequence space by column space: %d -> %d \n", sspace->capacity, N );

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
    sspace->seq_count = seq_count;
    sspace->seq_width = W;
}

template < class IntType >
__global__ void _delete_space( device_sequence_space< IntType > * sspace ) {
    if( sspace->sequences != NULL ) {
        delete sspace->sequences;
    }
}

#endif  // DEVICE_SEQUENCE_SPACE_KERNELS_HPP_
