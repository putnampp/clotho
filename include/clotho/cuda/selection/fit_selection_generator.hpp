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
#ifndef FIT_SELECTION_GENERATOR_HPP_
#define FIT_SELECTION_GENERATOR_HPP_

#include "clotho/cuda/curand_state_pool.hpp"

#include "clotho/cuda/data_spaces/basic_data_space.hpp"

#include "clotho/cuda/distributions/discrete_distribution.hpp"
//#include "clotho/cuda/selection/random_select_parents.hpp"
#include "clotho/cuda/recombination/recombine_parents.hpp"

#include "clotho/cuda/device_state_object.hpp"

template < typename IntType, typename RealType >
class FitSelectionGenerator : public clotho::utility::iStateObject {
public:
    typedef basic_data_space< IntType > event_space_type;
    
    typedef clotho::cuda::curand_state_pool             state_pool_type;

    typedef DiscreteDistribution< IntType, RealType >   discrete_dist_type;

    FitSelectionGenerator( boost::property_tree::ptree & config ) {
        state_pool_type::getInstance()->initialize( config );
        create_space( dParentIndices );
    }

    template < class PopulationSpaceType, class FitnessSpaceType >
    void generate( PopulationSpaceType * parent_pop, PopulationSpaceType * child_pop, FitnessSpaceType * fitness ) {
        m_discrete.initialize_table( fitness );
        CHECK_LAST_KERNEL_EXEC

        make_random_list<<< 1, 32 >>>( state_pool_type::getInstance()->get_device_states()
                                        , m_discrete.get_device_space()
                                        , parent_pop->sequences.get_device_space()
                                        , child_pop->sequences.get_device_space()
                                        , dParentIndices );
        CHECK_LAST_KERNEL_EXEC
    }

    template < class PopulationSpaceType, class FitnessSpaceType >
    void generate_and_recombine( PopulationSpaceType * parent_pop, PopulationSpaceType * child_pop, FitnessSpaceType * fitness ) {
        generate( parent_pop, child_pop, fitness );

        clotho::utility::algo_version< 2 > * v = NULL;
        recombine_parents_kernel<<< 10, 1024 >>>( parent_pop->sequences.get_device_space()
                                                , dParentIndices
                                                , child_pop->sequences.get_device_space() 
                                                , v);
        CHECK_LAST_KERNEL_EXEC
    }

    event_space_type * get_device_space() {
        return dParentIndices;
    }

    void get_state( boost::property_tree::ptree & state ) {
        boost::property_tree::ptree fit;
        m_discrete.get_state( fit );

        boost::property_tree::ptree evts;
        get_device_object_state( evts, dParentIndices );

        state.add_child( "distribution", fit );
        state.add_child( "events", evts );
    }

    virtual ~FitSelectionGenerator() {
        delete_space( dParentIndices );
    }

protected:
    event_space_type    * dParentIndices;
    discrete_dist_type  m_discrete;
};

#endif  // FIT_SELECTION_GENERATOR_HPP_
