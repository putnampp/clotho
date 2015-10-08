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

#include "clotho/cuda/data_spaces/sequence_space/device_sequence_space_def.hpp"
#include "clotho/cuda/data_spaces/sequence_space/device_sequence_space_helper.hpp"

#include "clotho/cuda/data_spaces/data_space_api.hpp"

#include "clotho/cuda/data_spaces/sequence_space/device_sequence_space_kernels.hpp"

template < class IntType, class SpaceType >
void resize_space( device_sequence_space< IntType > * seqs
                 , SpaceType * column_space
                 , unsigned int row_count ) {
    _resize_space<<< 1, 1 >>>( seqs, column_space, row_count );
}

template < class IntType >
void resize_space( device_sequence_space< IntType > * seqs
                    , unsigned int cols
                    , unsigned int rows ) {
    _resize_space<<< 1, 1 >>>( seqs, cols, rows );
}

#endif  // DEVICE_SEQUENCE_SPACE_HPP_
