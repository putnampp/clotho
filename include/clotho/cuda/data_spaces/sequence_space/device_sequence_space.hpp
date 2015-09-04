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
#ifndef DEVICE_SEQUENCE_SPACE_HPP_
#define DEVICE_SEQUENCE_SPACE_HPP_

#include "clotho/cuda/data_spaces/sequence_space/device_sequence_space_kernels.hpp"
#include "clotho/cuda/data_spaces/data_space_api.hpp"

/*
template < class SpaceType >
void create_sequence_space( SpaceType *& sspace ) {
    assert( cudaMalloc( (void **) &sspace, sizeof( SpaceType ) ) == cudaSuccess );

    assert( cudaMemset( sspace, 0, sizeof(SpaceType) ) == cudaSuccess );
}

template < class SSpaceType, class ASpaceType >
void resize_sequence_space( SSpaceType * sspace, ASpaceType * aspace, unsigned int seq_count ) {
    _resize_sequence_space<<< 1, 1 >>>( sspace, aspace, seq_count );
}

template < class SpaceType >
void delete_sequence_space( SpaceType * sspace ) {
    _delete_sequence_space<<< 1, 1 >>>( sspace );

    cudaFree( sspace );
}*/

#endif  // DEVICE_SEQUENCE_SPACE_HPP_
