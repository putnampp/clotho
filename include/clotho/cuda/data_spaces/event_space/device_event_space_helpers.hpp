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
#ifndef DEVICE_EVENT_SPACE_HELPERS_HPP_
#define DEVICE_EVENT_SPACE_HELPERS_HPP_

#include <cuda.h>
#include <ostream>

template < class IntType, class OrderTag >
void write_data( std::ostream & out, device_event_space< IntType, OrderTag > & lspace ) {
    out << "\n\"event_count\": []";
}

template < class IntType, class OrderTag >
void write_summary( std::ostream & out, device_event_space< IntType, OrderTag > & lspace ) {}

template < class IntType >
void write_summary( std::ostream & out, device_event_space< IntType, unit_ordered_tag< IntType > > & lspace ) {
    out << ",\n\"bin_summary\": [ " << lspace.bin_summary[0];
    for( unsigned int i = 1; i < unit_ordered_tag< IntType >::OBJECTS_PER_INT; ++i ) {
        out << ",\n" << lspace.bin_summary[i];
    }
    out << "]";
}

template < class IntType, class OrderTag >
void write_stats( std::ostream & out, device_event_space< IntType, OrderTag > & lspace ) {
    out << ",\n\"total\": " << lspace.total;
    out << ",\n\"size\": " << lspace.size;
    out << ",\n\"capacity\": " << lspace.capacity;
}

template < class IntType, class OrderTag >
std::ostream & operator<<( std::ostream & out, device_event_space< IntType, OrderTag > * space ){
    if( space == NULL ) return out;

    typedef device_event_space< IntType, OrderTag > space_type;

    space_type lspace;

    assert( cudaMemcpy( &lspace, space, sizeof( space_type ) , cudaMemcpyDeviceToHost ) == cudaSuccess );

    out << "{";
    write_data( out, lspace );
    write_summary( out, lspace );
    write_stats( out, lspace );
    out << "\n}";

    return out;
}

#endif  // DEVICE_EVENT_SPACE_HELPERS_HPP_
