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
#ifndef SEQUENCE_SPACE_HPP_
#define SEQUENCE_SPACE_HPP_

#include "clotho/cuda/data_spaces/sequence_space/device_sequence_space.hpp"

#include <algorithm>

template < class IntType >
class SequenceSpace {
public:
    typedef IntType int_type;
    typedef device_sequence_space< int_type > device_space_type;

    SequenceSpace() {
        initialize();
    }

    template < class ASpaceType >
    void resize_space( ASpaceType * aspace, unsigned int seq_count ) {
        resize_space( dSequences, aspace, seq_count );
    }

    virtual ~SequenceSpace() {
        delete_space( dSequences );
    }

protected:
    void initialize() {
        create_space( dSequences );
    }

    device_space_type * dSequences;
};

#endif  // SEQUENCE_SPACE_HPP_
