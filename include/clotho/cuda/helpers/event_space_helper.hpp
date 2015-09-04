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
#ifndef EVENT_SPACE_HELPER_HPP_
#define EVENT_SPACE_HELPER_HPP_

template < class OrderTag >
struct event_space_helper {

    static unsigned int get( unsigned int N ) {
        return N;
    }
};

#include "clotho/cuda/data_spaces/tags/unit_ordered_tag.hpp"
template < class IntType >
struct event_space_helper< unit_ordered_tag< IntType > > {
    static unsigned int get( unsigned int N ) {
        unsigned int n = N * unit_ordered_tag< IntType >::OBJECTS_PER_UNIT;
        return n;
    }
};

#endif  // EVENT_SPACE_HELPER_HPP_
