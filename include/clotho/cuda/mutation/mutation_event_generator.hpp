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

#include "clotho/cuda/curand_helper.hpp"

#include "clotho/cuda/distributions/poisson_distribution.hpp"

#include "clotho/cuda/mutation/event_space_helper.hpp"
#include "clotho/cuda/mutation/scaled_mean_helper.hpp"

template < class RealType, class IntType, class OrderTag >
class MutationEventGenerator {
public:
    typedef RealType  real_type;
    typedef device_event_space< IntType > space_type;

    typedef curandState_t                               state_type;
    typedef clotho::cuda::curand_helper< state_type >   helper_type;
    typedef typename helper_type::seed_type             seed_type;

    typedef poisson_cdf< real_type, 32 >                poisson_type;

    typedef scaled_mean_helper< real_type, OrderTag >                     scaled_mean_type;
    typedef event_space_helper< OrderTag >  space_helper_type;

    MutationEventGenerator( boost::property_tree::ptree & config ) :
        m_states( NULL )
        , m_seed( 0 )
        , m_mutation_rate( 0.001 )
    {
        parse_configuration( config );

        initialize();
    }

    void operator()( space_type * space, unsigned int N ) {
        unsigned int event_counts = space_helper_type::get( N );
        resize_space( space, event_counts );
        _simple_mutation_generator<<< 1, 32 >>>( m_states, space, dPoisCDF, N, (OrderTag * ) NULL);
    }

    virtual ~MutationEventGenerator() {
        helper_type::cleanup_states( m_states );
    }

protected:

    void initialize() {

        assert( cudaMalloc( (void **) &dPoisCDF, sizeof( poisson_type) ) == cudaSuccess );

        initialize_poisson( m_mutation_rate );

        helper_type::make_states( m_states, m_seed, 1, 32 );
    }

    void initialize_poisson( real_type mean ) {
        real_type smean = scaled_mean_type::get( mean );
        make_poisson_cdf_maxk32<<< 1, 32 >>>( dPoisCDF, smean );
    }

    real_type    get_scaled_mean( real_type m );

    unsigned int event_space_size( unsigned int N );

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

        if( lconfig.get_child_optional( "generator.seed" ) != boost::none ) {
            m_seed = lconfig.get< seed_type >( "generator.seed", m_seed );
        }

        if( m_seed == 0 ) {
            m_seed = clotho::utility::clock_type::now().time_since_epoch().count();
            lconfig.put("generator.seed", m_seed );
        }

        config.put_child( "mutation", lconfig );
    }

    state_type          * m_states;
    seed_type           m_seed;

    poisson_type        * dPoisCDF;
    
    real_type           m_mutation_rate;
};

#endif  // MUTATION_EVENT_GENERATOR_HPP_
