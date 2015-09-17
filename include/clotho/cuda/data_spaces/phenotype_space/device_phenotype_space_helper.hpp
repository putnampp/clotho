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
#ifndef DEVICE_PHENOTYPE_SPACE_HELPER_HPP_
#define DEVICE_PHENOTYPE_SPACE_HELPER_HPP_

#include "clotho/utility/state_object.hpp"

namespace clotho {
namespace utility {

template < class RealType >
//void get_state( boost::property_tree::ptree & state, const device_phenotype_space< RealType > & obj ) {
struct state_getter< device_phenotype_space< RealType > > {
    void operator()( boost::property_tree::ptree & state, const device_phenotype_space< RealType > & obj ) {
//        std::cerr << "Phenotype" << std::endl;
/*
        state.put( "size", obj.size );
        state.put( "capacity", obj.capacity );

        if( obj.data == NULL || obj.size == 0 ) return;

        typedef typename device_phenotype_space< RealType >::real_type real_type;

        real_type * ecount = new real_type[ obj.size ];

        copy_heap_data( ecount, obj.data, obj.size );

        boost::property_tree::ptree e;
        for( unsigned int i = 0; i < obj.size; ++i ) {
            clotho::utility::add_value_array( e, ecount[ i ] );
        }

        state.add_child( "data", e );

        delete ecount;
*/
        state_getter< basic_data_space< RealType > > sget;
        sget( state, (basic_data_space< RealType >) obj );
    }
};

}   // namespace utility
}   // namespace clotho
#endif  // DEVICE_PHENOTYPE_SPACE_HELPER_HPP_
