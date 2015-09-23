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
#ifndef DEVICE_FREE_SPACE_DEF_HPP_
#define DEVICE_FREE_SPACE_DEF_HPP_

template < class IntType, class OrderTag >
struct device_free_space {
    typedef IntType     int_type;
    typedef OrderTag    order_tag_type;

    static const unsigned int OBJECTS_PER_INT = (sizeof( int_type ) * 8);

    int_type        * lost_list;
    int_type        * fixed_list;
    int_type        * free_list;
    unsigned int    * free_map;

    unsigned int total, fixed_count;
    unsigned int size, capacity;
};

#endif  // DEVICE_FREE_SPACE_DEF_HPP_
