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
#ifndef POPULATION_MUTATION_GENERATOR_HPP_
#define POPULATION_MUTATION_GENERATOR_HPP_

#include <boost/property_tree/ptree.hpp>

#include "clotho/cuda/data_space/event_space/device_event_space.hpp"

#include "clotho/utility/timer.hpp"

#include <curand.h>
#include <curand_kernel.h>

#include "clotho/cuda/curand_helper.hpp"

template < class PopulationType >
class population_mutation_generator {
public:
    typedef PopulationType                              population_type;

    typedef typename population_type::sequence_space_type     sequence_space_type;

    typedef typename population_type::allele_space_type allele_space_type;

    typedef typename allele_space_type::real_type       real_type;
    typedef typename allele_space_type::order_tag_type  order_tag_type;
    typedef typename sequence_space_type::int_type      int_type;

    typedef device_event_space< int_type, order_tag_type >  event_space_type;

    typedef curandState_t                               state_type;
    typedef clotho::cuda::curand_helper< state_type >   helper_type;
    typedef typename helper_type::seed_type             seed_type;

    population_mutation_generator( boost::property_tree::ptree & config ) {
        parse_configuration( config );
        initialize();
    }

    void operator()( population_type * pop, event_space_type * mut ) {

    }

    virtual ~population_mutation_generator() {}

protected:
    void parse_configuration( boost::property_tree::ptree & config ) {
        boost::property_tree::ptree lconfig;

        if( config.get_child_optional( "mutation" ) != boost::none ) {
            lconfig = config.get_child( "mutation" );
        }

        if( lconfig.get_child_optional( "mutation_per_sequence" ) == boost::none ) {
            lconfig.put("mutation_per_sequence", m_mutation_rate );
        } else {
            m_mutation_rate = lconfig.get< real_type >( "mutation_per_sequence", m_mutation_rate );
        }

        if( lconfig.get_child_optional( "rng.seed" ) != boost::none ) {
            m_seed = lconfig.get< seed_type >( "rng.seed", m_seed );
        }

        if( m_seed == 0 ) {
            m_seed = clotho::utility::clock_type::now().time_since_epoch().count();
            lconfig.put("rng.seed", m_seed );
        }

        config.put_child( "mutation", lconfig );
    }

    state_type          * m_states;
    seed_type           m_seed;
};

#endif  // POPULATION_MUTATION_GENERATOR_HPP_
