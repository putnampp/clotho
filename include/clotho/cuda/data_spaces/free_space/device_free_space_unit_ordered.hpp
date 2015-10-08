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
#ifndef DEVICE_FREE_SPACE_UNIT_ORDERED_HPP_
#define DEVICE_FREE_SPACE_UNIT_ORDERED_HPP_

/*
#include "clotho/cuda/data_spaces/free_space/device_free_space_def.hpp"
#include "clotho/cuda/data_spaces/tags/unit_ordered_tag.hpp"

template < class IntType >
struct device_free_space< IntType, unit_ordered_tag< IntType > > {
    typedef IntType                     int_type;
    typedef unit_ordered_tag< IntType > ordered_tag_type;

    static const unsigned int OBJECTS_PER_INT = ordered_tag_type::OBJECTS_PER_UNIT;

    int_type    * free_list;
    int_type    * free_map;

    unsigned int free_count[ ordered_tag_type::OBJECTS_PER_UNIT ];

    unsigned int total;
    unsigned int size, capacity;
};*/

#endif  // DEVICE_FREE_SPACE_UNIT_ORDERED_HPP_
