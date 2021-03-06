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
#ifndef DEVICE_FREE_SPACE_HELPER_HPP_
#define DEVICE_FREE_SPACE_HELPER_HPP_

#include <ostream>

#include "clotho/cuda/data_spaces/data_space_helper.hpp"

template < class IntType, class OrderTag >
void dump_free_list( std::ostream & out, const device_free_space< IntType, OrderTag > & lspace ) {
    if( lspace.free_list == NULL || lspace.size == 0 ) {
        out << "[]";
        return;
    }

    typedef device_free_space< IntType, OrderTag >  space_type;
    typedef typename space_type::int_type           int_type;

    unsigned int N = lspace.size / space_type::OBJECTS_PER_INT;

    int_type    * _list = new int_type[ N ];

    copy_heap_data( _list, lspace.free_list, N );

    out << "[0x" << std::hex << std::setw(sizeof(int_type) * 2) << std::setfill('0') << _list[0];
    for( unsigned int i = 1; i < N; ++i ) {
        out << ",0x" << std::hex << std::setw(sizeof(int_type) * 2) << std::setfill('0') << _list[i];
    }

    out << "]";
    out << std::dec;

    delete _list;
}

template < class IntType, class OrderTag >
void dump_free_map( std::ostream & out, const device_free_space< IntType, OrderTag > & lspace ) {
    if( lspace.free_map == NULL || lspace.size == 0) {
        out << "[]";
        return;
    }

    typedef unsigned int int_type;

    unsigned int N = lspace.size;

    int_type    * _list = new int_type[ N ];

    copy_heap_data( _list, lspace.free_map, N );

    out << "[" << _list[0];
    for( unsigned int i = 1; i < N; ++i ) {
        out << "," << _list[i];
    }

    out << "]";

    delete _list;
}

template < class IntType, class OrderTag >
std::ostream & operator<<( std::ostream & out, device_free_space< IntType, OrderTag > * fspace ) {
    if( fspace == NULL ) {
        out << "{}";
        return out;
    }

    device_free_space< IntType, OrderTag > loc;

    assert( cudaMemcpy( &loc, fspace, sizeof( device_free_space< IntType, OrderTag > ), cudaMemcpyDeviceToHost ) == cudaSuccess );

    out << "{";
    out << "\n\"free_list\": ";
#ifdef DUMP_HEAP_VARIABLES
    dump_free_list( out, loc );
#else
    out << "[]";
#endif  // DUMP_HEAP_VARIABLES

    out << ",\n\"free_map\": ";
#ifdef DUMP_HEAP_VARIABLES
    dump_free_map( out, loc );
#else
    out << "[]";
#endif  // DUMP_HEAP_VARIABLES

    out << ",\n\"total\": " << loc.total;
    out << ",\n\"size\": " << loc.size;
    out << ",\n\"capacity\": " << loc.capacity;
    out << "\n}";

    return out;
}

#include "clotho/utility/state_object.hpp"

namespace clotho {
namespace utility {

template < class IntType, class OrderTag >
//void get_state( boost::property_tree::ptree & state, const device_free_space< IntType, OrderTag > & obj ) {
struct state_getter< device_free_space< IntType, OrderTag > > {

    void operator()( boost::property_tree::ptree & state, const device_free_space< IntType, OrderTag > & obj ) {
//        state.put( "total", obj.total);
//
        typedef device_free_space< IntType, OrderTag > space_type;
        typedef typename space_type::int_type  int_type;

        state_getter < device_event_space< IntType, OrderTag > > base_getter;
        base_getter( state, obj );

        state.put( "fixed_count", obj.fixed_count );
        state.put( "size", obj.size );
        state.put( "capacity", obj.capacity );

        if( obj.free_list == NULL ) {
            state.put( "free_list", "" );
            state.put( "free_map", "" );
            state.put( "fixed_list", "" );
            state.put( "lost_list", "" );
            return;
        }

        unsigned int nInts = obj.size / space_type::OBJECTS_PER_INT;
        int_type * flist = new int_type[ nInts ];

        copy_heap_data( flist, obj.free_list, nInts );

        boost::property_tree::ptree fl, xl, ll;
        for( unsigned int i = 0; i < nInts; ++i ) {
            std::ostringstream oss;
            oss << "0x" << std::hex << std::setw( sizeof( int_type) * 2 ) << std::setfill( '0' ) << flist[i];

            clotho::utility::add_value_array( fl, oss.str() );
        }

        copy_heap_data( flist, obj.fixed_list, nInts );

        for( unsigned int i = 0; i < nInts; ++i ) {
            std::ostringstream oss;
            oss << "0x" << std::hex << std::setw( sizeof( int_type) * 2 ) << std::setfill( '0' ) << flist[i];

            clotho::utility::add_value_array( xl, oss.str() );
        }

        copy_heap_data( flist, obj.lost_list, nInts );

        for( unsigned int i = 0; i < nInts; ++i ) {
            std::ostringstream oss;
            oss << "0x" << std::hex << std::setw( sizeof( int_type) * 2 ) << std::setfill( '0' ) << flist[i];

            clotho::utility::add_value_array( ll, oss.str() );
        }

        unsigned int * fmap = new unsigned int[ obj.size ];

        copy_heap_data( fmap, obj.free_map, obj.size );

        boost::property_tree::ptree fm;
        for( unsigned int i = 0; i < obj.size; ++i ) {
            clotho::utility::add_value_array( fm, fmap[i] );
        }

        state.add_child( "free_list", fl );
        state.add_child( "free_map", fm );
        state.add_child( "lost_list", ll );
        state.add_child( "fixed_list", xl );

        delete flist;
        delete fmap;
    }
};

}   // namespace utility
}   // namespace clotho

#endif  // DEVICE_FREE_SPACE_HELPER_HPP_
