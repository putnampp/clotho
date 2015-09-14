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

#include "clotho/cuda/data_spaces/data_space_helper.hpp"
#include "clotho/utility/log_helper.hpp"

template < class IntType, class OrderTag >
void write_data( std::ostream & out, device_event_space< IntType, OrderTag > & lspace ) {

    if( lspace.size == 0 ) {
        out << "[]";
        return;
    }

    typedef device_event_space< IntType, OrderTag > space_type;
    typedef typename space_type::int_type           event_int_type;

    event_int_type * events = new event_int_type[ lspace.size ];

    copy_heap_data( events, lspace.event_count, lspace.size );

    out << "[" << events[0];
    for( int i = 1; i < lspace.size; ++i ) {
        out << ",\n" << events[i];
    }

    out << "]";
//    cudaFree( dEvt );

    delete events;
}

template < class OutputType, class IntType, class OrderTag >
void write_summary( OutputType & out, const device_event_space< IntType, OrderTag > & lspace ) {}

template < class IntType >
void write_summary( std::ostream & out, const device_event_space< IntType, unit_ordered_tag< IntType > > & lspace ) {
    out << ",\n\"bin_summary\": [ " << lspace.bin_summary[0];
    for( unsigned int i = 1; i < unit_ordered_tag< IntType >::OBJECTS_PER_UNIT; ++i ) {
        out << ",\n" << lspace.bin_summary[i];
    }
    out << "]";
}

template < class IntType >
void write_summary( boost::property_tree::ptree & out, const device_event_space< IntType, unit_ordered_tag< IntType > > & lspace ) {
    typedef device_event_space< IntType, unit_ordered_tag< IntType > > space_type;
    
    boost::property_tree::ptree s;
    for( unsigned int i = 0; i < space_type::order_tag_type::OBJECTS_PER_UNIT; ++i ) {
        clotho::utility::add_value_array( s, lspace.bin_summary[i] );
    }

    out.add_child( "bin_summary", s );
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
    out << "\"event_count\": ";
#ifdef DUMP_HEAP_VARIABLES
    write_data( out, lspace );
#else
    out << "[]";
#endif  // DUMP_HEAP_VARIABLES
    write_summary( out, lspace );
    write_stats( out, lspace );
    out << "\n}";

    return out;
}

#include "clotho/utility/state_object.hpp"

namespace clotho {
namespace utility {

//template < class IntType, class OrderTag >
//void get_state( boost::property_tree::ptree & state, const device_event_space< IntType, OrderTag > & obj ) {

template < class IntType, class OrderTag >
struct state_getter< device_event_space< IntType, OrderTag > > {
    void operator()( boost::property_tree::ptree & state, const device_event_space< IntType, OrderTag > & obj ) {
        state.put( "total", obj.total );
        state.put( "size", obj.size );
        state.put( "capacity", obj.capacity );

        typedef typename device_event_space< IntType, OrderTag >::int_type int_type;

        int_type * ecount = new int_type[ obj.size ];

        copy_heap_data( ecount, obj.event_count, obj.size );

        boost::property_tree::ptree e;
        for( unsigned int i = 0; i < obj.size; ++i ) {
            clotho::utility::add_value_array( e, ecount[ i ] );
        }

        state.add_child( "event_count", e );

        write_summary( state, obj );

        delete ecount;
    }
};

}   // namespace utility
}   // namespace clotho

#endif  // DEVICE_EVENT_SPACE_HELPERS_HPP_
