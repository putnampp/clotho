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
#ifndef DEVICE_ALLELE_SPACE_HELPER_HPP_
#define DEVICE_ALLELE_SPACE_HELPER_HPP_

#include "clotho/cuda/data_spaces/data_space_helper.hpp"

//template < class DataType >
//__global__ void copy_heap( DataType * d_loc, DataType * d_heap, unsigned int N ) {
//    memcpy( d_loc, d_heap, N * sizeof( DataType ) );
//}

template < class RealType/*, class IntType, class OrderTag*/ >
void dump_locations( std::ostream & out, const device_allele_space< RealType/*, IntType, OrderTag*/ > & rhs ) {
    if( rhs.locations == NULL || rhs.size == 0 ) {
        out << "[]";
        return;
    }

    typedef device_allele_space< RealType/*, IntType, OrderTag*/ >  space_type;
    typedef typename space_type::real_type                      real_type;

    unsigned int N = rhs.size;
    real_type * locations = new real_type[ N ];

/*    real_type * dLoc;
    assert( cudaMalloc( (void **) &dLoc, N * sizeof( real_type ) ) == cudaSuccess );

    copy_heap<<< 1, 1 >>>( dLoc, rhs.locations, N );
    cudaDeviceSynchronize();

    assert( cudaMemcpy( locations, dLoc, N * sizeof( real_type ), cudaMemcpyDeviceToHost ) == cudaSuccess );
*/
    copy_heap_data( locations, rhs.locations, N );

    out << "[" << locations[0];

    for( int i = 1; i < N; ++i ) {
        out << ",\n" << locations[i];
    }
    out << "]";

//    cudaFree( dLoc );
    delete locations;
}

/*
template < class RealType, class IntType, class OrderTag >
void dump_free_list( std::ostream & out, const device_allele_space< RealType, IntType, OrderTag > & rhs ) {
    if( rhs.free_list == NULL || rhs.size == 0 ) {
        out << "[]";
        return;
    }

    typedef typename device_allele_space< RealType, IntType, OrderTag >::int_type int_type;
    unsigned int N = compute_free_list_size< IntType >( rhs.size );
    int_type * free_list = new int_type[ N ];

    int_type * d_free;
    assert( cudaMalloc( (void **) &d_free, N * sizeof( int_type ) ) == cudaSuccess );

    copy_heap<<< 1, 1>>>( d_free, rhs.free_list, N );
    assert( cudaMemcpy( free_list, d_free, N * sizeof( int_type ), cudaMemcpyDeviceToHost ) == cudaSuccess );

    out << std::hex;
    out << "[0x" << std::setw(sizeof(int_type) * 2) << std::setfill('0') << free_list[0];

    for( int i = 1; i < N; ++i ) {
        out << ",0x" << std::setw(sizeof(int_type) * 2) << std::setfill('0') << free_list[i];
    }
    out << "]";

    out << std::dec;

    cudaFree( d_free );
    delete free_list;
}

template < class RealType, class IntType, class OrderTag >
void dump_free_count( std::ostream & out, const device_allele_space< RealType, IntType, OrderTag > & rhs ) {
    out << rhs.free_count;
}

template < class RealType, class IntType >
void dump_free_count ( std::ostream & out, const device_allele_space< RealType, IntType, unit_ordered_tag< IntType > > & rhs ) {
    out << "[" << rhs.free_count[0];

    for( int i = 1; i < unit_ordered_tag< IntType >::OBJECTS_PER_UNIT; ++i ) {
        out << ",\n" << rhs.free_count[i];
    }
    out << "]";
}*/

template < class RealType/*, class IntType, class OrderTag*/ >
std::ostream & operator<<( std::ostream & out, const device_allele_space< RealType/*, IntType, OrderTag*/ > & rhs ) {
    out << "{";

    out << "\n\"locations\": "; // << std::hex << rhs.locations;
#ifdef DUMP_HEAP_VARIABLES
    dump_locations( out, rhs );
#else
    out << "[]";
#endif  // DUMP_HEAP_VARIABLES

//    out << ",\n\"free_list\": ";
//#ifdef DUMP_HEAP_VARIABLES
//    dump_free_list( out, rhs );
//#else
//    out << "[]";
//#endif  // DUMP_HEAP_VARIABLES
//    out << ",\n\"free_count\": ";
//    dump_free_count( out, rhs );

    out << ",\n\"size\": " << rhs.size;
    out << ",\n\"capacity\": " << rhs.capacity;
    out << "}";
    return out;
}

#include "clotho/utility/state_object.hpp"
#include "clotho/utility/log_helper.hpp"

namespace clotho {
namespace utility {

//template < class RealType > 
//void get_state( boost::property_tree::ptree & state, const device_allele_space< RealType > & obj ) {

template < class RealType >
struct state_getter< device_allele_space< RealType > > {

    void operator()( boost::property_tree::ptree & state, const device_allele_space< RealType > & obj ) {
//        std::cerr << "Basic allele" << std::endl;

        //state.put( "size", obj.size );
        state.put( "capacity", obj.capacity );

        unsigned int N = obj.capacity;
        if( N == 0 ) return;

        typedef device_allele_space< RealType > space_type;
        typename space_type::real_type * loc = new typename space_type::real_type[ N ];

        copy_heap_data( loc, obj.locations, N );

        boost::property_tree::ptree l;
        for( unsigned int i = 0; i < N; ++i ) {
            clotho::utility::add_value_array( l, loc[i] );
        }

        state.add_child( "locations", l );
        delete loc;
    }
};

//template < class RealType >
//void get_state( boost::property_tree::ptree & state, const device_weighted_allele_space< RealType > & obj ) {
template < class RealType >
struct state_getter< device_weighted_allele_space < RealType > > {

    void operator()( boost::property_tree::ptree & state, const device_weighted_allele_space< RealType > & obj ) {
        state_getter< device_allele_space< RealType > > base_getter;
        base_getter( state, (device_allele_space< RealType >) obj );

        //std::cerr<< "Weighted allele" << std::endl;

        unsigned int N = obj.capacity;
        if( N == 0 ) return;

        typedef typename device_weighted_allele_space< RealType >::real_type real_type;
        real_type   * weights = new real_type[ N ];

        copy_heap_data( weights, obj.weights, N );

        boost::property_tree::ptree w;
        for( unsigned int i = 0; i < N; ++i ) {
            clotho::utility::add_value_array( w, weights[i] );
        }

        state.add_child( "weights", w );

        delete weights;
    }
};

}   // namespace utility
}   // namespace clotho

#endif  // DEVICE_ALLELE_SPACE_HELPER_HPP_
