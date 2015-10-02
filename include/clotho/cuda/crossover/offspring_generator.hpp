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
#ifndef OFFSPRING_GENERATOR_HPP_
#define OFFSPRING_GENERATOR_HPP_

#include <boost/property_tree/ptree.hpp>

#include "clotho/cuda/curand_state_pool.hpp"

#include "clotho/cuda/crossover/select_and_crossover_kernels.hpp"

#include "clotho/cuda/distributions/poisson_distribution.hpp"
#include "clotho/recombination/recombination_rate_parameter.hpp"

template < class PopulationType >
class OffspringGenerator {
public:

    typedef PopulationType                                  population_type;
    typedef typename population_type::real_type             real_type;
    typedef typename population_type::int_type              int_type;
    typedef typename population_type::order_tag_type        order_tag_type;

    typedef clotho::cuda::curand_state_pool                 state_pool_type;

    typedef recombination_rate_parameter< real_type > recombination_rate_type;
    typedef poisson_cdf< real_type, 32 >                    poisson_type;

    OffspringGenerator( boost::property_tree::ptree & config ) :
        m_recomb_rate( config )
    {
        parse_configuration( config );
        initialize();
    }

    void operator()( population_type * parent, selection_type * fit, population_type * offspring ) {

        unsigned int bcount = state_pool_type::getInstance()->get_max_blocks();
        unsigned int tcount = state_pool_type::getInstance()->get_max_threads();

        select_and_crossover_kernel<<< bcount, tcount >>>( state_pool_type::getInstance()->get_device_states(), pop->sequences.get_device_space(), 
    }

    virtual ~OffspringGenerator() {
        cudaFree( dPoisCDF );
    }

protected:

    void parse_configuration( boost::property_tree::ptree & config ) {
        state_pool_type::getInstance()->initialize( config );
    }

    void initialize() {
        assert( cudaMalloc( (void **) &dPoisCDF, sizeof( poisson_type) ) == cudaSuccess );

        initialize_poisson( m_recombination_rate.m_rho, (order_tag_type *) NULL );
    }
    poisson_type            * dPoisCDF;
    recombination_rate_type m_recomb_rate;
};

#endif  // OFFSPRING_GENERATOR_HPP_
