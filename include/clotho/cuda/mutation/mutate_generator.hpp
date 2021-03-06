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
#ifndef MUTATE_GENERATOR_HPP_
#define MUTATE_GENERATOR_HPP_

#include <boost/property_tree/ptree.hpp>

#include "clotho/cuda/data_spaces/event_space/device_event_space.hpp"
#include "clotho/cuda/mutation/simple_mutation_generator.hpp"
#include "clotho/utility/timer.hpp"

#include <curand.h>
#include <curand_kernel.h>

#include "clotho/cuda/curand_state_pool.hpp"

#include "clotho/cuda/distributions/poisson_distribution.hpp"

#include "clotho/cuda/mutation/generate_mutation_kernel_api.hpp"
#include "clotho/cuda/mutation/generate_mutation_kernel_impl.hpp"

#include "clotho/cuda/mutation/scatter_unordered_impl.hpp"

#include "clotho/cuda/device_state_object.hpp"
#include "clotho/mutation/mutation_rate_parameter.hpp"

template < class RealType, class IntType, class OrderTag >
class MutateGenerator : public clotho::utility::iStateObject {
public:
    typedef RealType  real_type;
    typedef device_event_space< IntType, OrderTag > space_type;

    typedef clotho::cuda::curand_state_pool             state_pool_type;

    typedef poisson_cdf< real_type, 32 >                poisson_type;

    //typedef scaled_mean_helper< real_type, OrderTag >                     scaled_mean_type;
    //typedef event_space_helper< OrderTag >  space_helper_type;

    MutateGenerator( boost::property_tree::ptree & config ) :
        m_mutation_rate( config )
    {
        parse_configuration( config );

        initialize();
    }

    void generate( space_type * space, unsigned int N ) {
        make_event_distribution_kernel<<< 1, 32 >>>( state_pool_type::getInstance()->get_device_space(), space, m_mutation_rate.m_mu, N );
    }

    template < class PopulationType >
    void scatter( PopulationType * pop, space_type * events, unsigned int N ) {
        const unsigned int MAX_BLOCKS = 40000;  // arbitrary limitation (think 65535 is max for any single grid dimension)
        unsigned int offset = 0;

        while( offset < N ) {
            unsigned int bcount = N - offset;
            bcount = (( bcount > MAX_BLOCKS ) ? MAX_BLOCKS : bcount);

            _scatter_mutation_single_thread<<< bcount, 1 >>>( pop->free_space, events, pop->sequences.get_device_space(), offset );
            offset += bcount;
        }

        _generate_mutation_kernel<<< 1, 32 >>>( state_pool_type::getInstance()->get_device_states(), pop->free_space, events, pop->alleles.get_device_space() );
    }

    void get_state( boost::property_tree::ptree & state ) {
    }

    virtual ~MutateGenerator() {
        cudaFree( dPoisCDF );
    }

protected:

    void initialize() {

        assert( cudaMalloc( (void **) &dPoisCDF, sizeof( poisson_type) ) == cudaSuccess );

        initialize_poisson( m_mutation_rate.m_mu );
    }

    void initialize_poisson( real_type mean ) {
        real_type smean = scaled_mean_type::get( mean );
        make_poisson_cdf_maxk32<<< 1, 32 >>>( dPoisCDF, smean );
    }

    void parse_configuration( boost::property_tree::ptree & config ) {
        state_pool_type::getInstance()->initialize(config);
    }

    poisson_type        * dPoisCDF;
    
    mutation_rate_parameter< real_type > m_mutation_rate;
};

#endif  // MUTATE_GENERATOR_HPP_
