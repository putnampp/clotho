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
#ifndef DEVICE_EVENT_SPACE_DEF_HPP_
#define DEVICE_EVENT_SPACE_DEF_HPP_

template < class IntType >
struct device_event_space {
    typedef IntType     int_type;

    typedef int_type *  pointer;

    pointer         event_count;

    unsigned int    total;
    unsigned int    size, capacity;
};

#include <cuda.h>
#include <ostream>

template < class IntType >
std::ostream & operator<<( std::ostream & out, device_event_space< IntType > * space ){
    if( space == NULL ) return out;

    typedef device_event_space< IntType > space_type;

    space_type lspace;

    assert( cudaMemcpy( &lspace, space, sizeof( space_type ) , cudaMemcpyDeviceToHost ) == cudaSuccess );

    out << "{";
    out << "\n\"total\":" << lspace.total;
    out << "\n}";

    return out;
}

#endif  // DEVICE_EVENT_SPACE_DEF_HPP_
