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
#ifndef MUTATION_EVENT_GENERATOR_HPP_
#define MUTATION_EVENT_GENERATOR_HPP_

#include <boost/property_tree/ptree.hpp>

#include "clotho/cuda/data_spaces/event_space/device_event_space.hpp"
#include "clotho/cuda/mutation/simple_mutation_generator.hpp"
#include "clotho/utility/timer.hpp"

#include <curand.h>
#include <curand_kernel.h>

#include "clotho/cuda/curand_state_pool.hpp"

#include "clotho/cuda/mutation/generate_mutation_kernel_api.hpp"
#include "clotho/cuda/mutation/generate_mutation_kernel_impl.hpp"

#include "clotho/cuda/mutation/scatter_unordered_impl.hpp"

#include "clotho/cuda/device_state_object.hpp"
#include "clotho/mutation/mutation_rate_parameter.hpp"

template < class RealType, class IntType, class OrderTag >
class MutationEventGenerator : public clotho::utility::iStateObject {
public:
    typedef RealType  real_type;
    typedef device_event_space< IntType, OrderTag > space_type;

    typedef clotho::cuda::curand_state_pool             state_pool_type;

    MutationEventGenerator( boost::property_tree::ptree & config ) :
        m_mutation_rate( config )
    {
        parse_configuration( config );
    }

/**
 * N - the number of sequences in the population
 */
    void generate( space_type * space, unsigned int N ) {
        make_event_distribution_kernel<<< 1, 32 >>>( state_pool_type::getInstance()->get_device_states(), space, m_mutation_rate.m_mu, N );
    }

    template < class PopulationType >
    void scatter( PopulationType * pop, space_type * events, unsigned int N ) {
        // scatter combines scattering and generation of allele
        //
        dim3 blocks( state_pool_type::getInstance()->get_total_states(), 1, 1), threads( 1,1,1);

        _scatter_mutations<<< blocks, threads >>>( state_pool_type::getInstance()->get_device_states()
                                            , pop->free_space
                                            , events
                                            , pop->sequences.get_device_space() 
                                            , pop->alleles.get_device_space() );
    }

    void get_state( boost::property_tree::ptree & state ) {
    }

    virtual ~MutationEventGenerator() { }

protected:

    void parse_configuration( boost::property_tree::ptree & config ) {
        state_pool_type::getInstance()->initialize(config);
    }

    mutation_rate_parameter< real_type > m_mutation_rate;
};

#endif  // MUTATION_EVENT_GENERATOR_HPP_
