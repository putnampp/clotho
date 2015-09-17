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
#ifndef BASIC_DATA_SPACE_HELPERS_HPP_
#define BASIC_DATA_SPACE_HELPERS_HPP_


#include "clotho/utility/state_object.hpp"

namespace clotho {
namespace utility {

template < class ValueType >
struct state_getter< basic_data_space< ValueType > > {
    void operator()( boost::property_tree::ptree & state, const basic_data_space< ValueType > & obj ) {
        state.put( "size", obj.size );
        state.put( "capacity", obj.capacity );

        if( obj.data == NULL || obj.size == 0 ) {
            state.put( "data", "");
            return;
        }

        typedef typename basic_data_space< ValueType >::value_type value_type;
        
        value_type  * dat = new value_type[ obj.size ];

        copy_heap_data( dat, obj.data, obj.size );

        boost::property_tree::ptree d;
        for( unsigned int i = 0; i < obj.size; ++i ) {
            clotho::utility::add_value_array( d, dat[i] );
        }
        state.add_child( "data", d );

        delete dat;
    }
};

}   // namespace utility
}   // namespace clotho
#endif  // BASIC_DATA_SPACE_HELPERS_HPP_
