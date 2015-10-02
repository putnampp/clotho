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
#ifndef CROSSOVER_GENERATOR_HPP_
#define CROSSOVER_GENERATOR_HPP_

#include <boost/property_tree/ptree.hpp>

#include "clotho/cuda/crossover/crossover_kernels.hpp"

#include "clotho/cuda/curand_state_pool.hpp"

#include "clotho/cuda/distributions/poisson_distribution.hpp"
#include "clotho/recombination/recombination_rate_parameter.hpp"

template < class PopulationType >
class CrossoverGenerator {
public:
    typedef PopulationType                                                      population_type;
    typedef typename population_type::allele_space_type::device_space_type      allele_space_type;

    typedef typename population_type::real_type             real_type;
    typedef typename population_type::int_type              int_type;
    typedef typename population_type::order_tag_type        order_tag_type;

    typedef clotho::cuda::curand_state_pool                 state_pool_type;

    typedef recombination_rate_parameter< real_type >       recombination_rate_type;

    typedef poisson_cdf< real_type, 32 >                    poisson_type;
    typedef scaled_mean_helper< real_type, order_tag_type > scaled_mean_type;

    CrossoverGenerator( boost::property_tree::ptree & config ) :
        //m_recombination_rate( 0.001 )
        m_recombination_rate( config )
    {
        parse_configuration( config );

        initialize();
    }

    void operator()( population_type * pop ) {
        unsigned int bcount = state_pool_type::getInstance()->get_max_blocks();
        unsigned int tcount = state_pool_type::getInstance()->get_max_threads();
        crossover_kernel<<< bcount, tcount >>>( state_pool_type::getInstance()->get_device_states(), (device_allele_space< real_type > * ) pop->alleles.get_device_space(), pop->free_space, dPoisCDF, pop->sequences.get_device_space() );
    }

    virtual ~CrossoverGenerator() {
        cudaFree( dPoisCDF );
    }

protected:

    void parse_configuration( boost::property_tree::ptree & config ) {
//        boost::property_tree::ptree lconfig;
//
//        if( config.get_child_optional( "recombination" ) != boost::none ) {
//            lconfig = config.get_child( "recombination" );
//        }
//
//        if( lconfig.get_child_optional( "crossover_per_sequence" ) == boost::none ) {
//            lconfig.put("crossover_per_sequence", m_recombination_rate );
//        } else {
//            m_recombination_rate = lconfig.get< real_type >( "crossover_per_sequence", m_recombination_rate );
//        }
        

        state_pool_type::getInstance()->initialize( config );

//        config.put_child( "recombination", lconfig );
    }

    void initialize() {

        assert( cudaMalloc( (void **) &dPoisCDF, sizeof( poisson_type) ) == cudaSuccess );

        initialize_poisson( m_recombination_rate.m_rho, (order_tag_type *) NULL );
    }

    void initialize_poisson( real_type mean, unordered_tag * p ) {
        real_type lambda = (mean / (real_type) allele_space_type::ALIGNMENT_SIZE);

        make_poisson_cdf_maxk32<<< 1, 32 >>>( dPoisCDF, lambda );
    }

    void initialize_poisson( real_type mean, unit_ordered_tag< int_type > * p ) {
        real_type lambda = (mean / 32.0);

        make_poisson_cdf_maxk32<<< 1, 32 >>>( dPoisCDF, lambda );
    }

    poisson_type * dPoisCDF;
    recombination_rate_type  m_recombination_rate;
};

#endif  // CROSSOVER_GENERATOR_HPP_
