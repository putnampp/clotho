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
#include "clotho/recombination/sequence_bias_parameter.hpp"

#include "clotho/cuda/distributions/bernoulli_distribution.hpp"

#include "clotho/cuda/device_state_object.hpp"

#include "clotho/cuda/selection/fit_selection_generator_def.hpp"

struct cuda_discrete_distribution {};

template < class IntType, class RealType >
class FitSelectionGenerator< IntType, RealType, cuda_discrete_distribution > : public clotho::utility::iStateObject {
public:
    typedef basic_data_space< IntType > event_space_type;
    
    typedef clotho::cuda::curand_state_pool             state_pool_type;

    typedef DiscreteDistribution< IntType, RealType >   discrete_dist_type;
    typedef sequence_bias_parameter< RealType >         sequence_bias_type;

    FitSelectionGenerator( boost::property_tree::ptree & config ) :
        dParentIndices(NULL)
        , m_seq_bias( config )
    {
        state_pool_type::getInstance()->initialize( config );
        create_space( dParentIndices );
    }

    template < class PopulationSpaceType, class FitnessSpaceType >
    void generate( PopulationSpaceType * parent_pop, PopulationSpaceType * child_pop, FitnessSpaceType * fitness, unsigned int N ) {
        m_discrete.initialize_table( fitness );
        CHECK_LAST_KERNEL_EXEC

        make_random_list<<< 1, 32 >>>( state_pool_type::getInstance()->get_device_states()
                                        , m_discrete.get_device_space()
                                        , parent_pop->sequences.get_device_space()
                                        , child_pop->sequences.get_device_space()
                                        , dParentIndices );
        CHECK_LAST_KERNEL_EXEC

        int shift = 1;
        inline_bernoulli_linear_shift_kernel<<< 1, 32 >>>( state_pool_type::getInstance()->get_device_states(), dParentIndices, m_seq_bias.m_bias, shift );
    }

    template < class PopulationSpaceType, class FitnessSpaceType >
    void generate_and_recombine( PopulationSpaceType * parent_pop, PopulationSpaceType * child_pop, FitnessSpaceType * fitness, unsigned int N ) {
        generate( parent_pop, child_pop, fitness, N );

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

    sequence_bias_type  m_seq_bias;
    discrete_dist_type  m_discrete;
};

#endif  // FIT_SELECTION_GENERATOR_HPP_
