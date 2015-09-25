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
#ifndef SELECTION_EVENT_GENERATOR_HPP_
#define SELECTION_EVENT_GENERATOR_HPP_

#include <boost/property_tree/ptree.hpp>

#include "clotho/cuda/curand_state_pool.hpp"

#include "clotho/cuda/data_spaces/event_space/device_event_space.hpp"
#include "clotho/cuda/data_spaces/tags/no_order_tag.hpp"

#include "clotho/cuda/selection/random_select_parents.hpp"
#include "clotho/cuda/recombination/recombine_parents.hpp"

#include "clotho/cuda/device_state_object.hpp"

class SelectionEventGenerator : public clotho::utility::iStateObject {
public:
    typedef device_event_space< unsigned int, no_order_tag > event_space_type;

    typedef clotho::cuda::curand_state_pool             state_pool_type;


    SelectionEventGenerator( boost::property_tree::ptree & config ) {
        parse_configuration( config );
        initialize();
    }

    template < class PopulationType >
    void generate( PopulationType * parent_pop, PopulationType * child_pop ) {
        typename state_pool_type::state_type * s = state_pool_type::getInstance()->get_device_states();

        random_select_parents_kernel<<< 1, 32 >>>( s
                                    , parent_pop->sequences.get_device_space()
                                    , child_pop->sequences.get_device_space()
                                    , dEvents );

        recombine_parents_kernel<<< 200, 32 >>>( parent_pop->sequences.get_device_space()
                                                , dEvents
                                                , child_pop->sequences.get_device_space() );
    }

    void get_state( boost::property_tree::ptree & state ) {
        get_device_object_state( state, dEvents );
    }

    virtual ~SelectionEventGenerator() {
        delete_space( dEvents );
    }
protected:

    void parse_configuration( boost::property_tree::ptree & config ) {
        state_pool_type::getInstance()->initialize( config );
    }

    void initialize() {
        create_space( dEvents );
    }

    event_space_type    * dEvents;   
};

#endif  // SELECTION_EVENT_GENERATOR_HPP_

