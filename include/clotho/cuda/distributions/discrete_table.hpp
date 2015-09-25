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
#ifndef DISCRETE_TABLE_HPP_
#define DISCRETE_TABLE_HPP_

template < class IntType, class RealType >
struct discrete_table {
    typedef RealType    real_type;
    typedef IntType     int_type;

    real_type   * threshold;
    int_type    * alternative;

    unsigned int size, capacity;
};

#include "clotho/cuda/distributions/discrete_table_kernels.hpp"

#include "clotho/utility/state_object.hpp"

namespace clotho {
namespace utility {

template < class IntType, class RealType >
struct state_getter< discrete_table< IntType, RealType > > {

    void operator()( boost::property_tree::ptree & state, const discrete_table< IntType, RealType > & obj ) {
        typedef discrete_table< IntType, RealType > table_type;
        typedef typename table_type::real_type      real_type;
        typedef typename table_type::int_type       int_type;

        state.put( "size", obj.size );
        state.put( "capacity", obj.capacity );

        if( obj.size == 0 || obj.threshold == NULL ) return;

        real_type * thresh = new real_type[ obj.size ];
        int_type * alts = new int_type[ obj.size ];

        copy_heap_data( thresh, obj.threshold, obj.size );
        copy_heap_data( alts, obj.alternative, obj.size );

        boost::property_tree::ptree t, a;
        for( unsigned int i = 0; i < obj.size; ++i ) {
            clotho::utility::add_value_array( t, thresh[i] );
        }
        state.add_child( "threshold", t );

        for( unsigned int i = 0; i < obj.size; ++i ) {
            clotho::utility::add_value_array( a, alts[i] );
        }
        state.add_child( "alternative", a );

        delete thresh;
        delete alts;
    }
};

}   // namespace utility
}   // namespace clotho

#endif  // DISCRETE_TABLE_HPP_
