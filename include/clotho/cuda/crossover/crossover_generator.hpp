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

//#include <curand.h>
//#include <curand_kernel.h>

//#include "clotho/cuda/curand_helper.hpp"
#include "clotho/cuda/curand_state_pool.hpp"

#include "clotho/cuda/distributions/poisson_distribution.hpp"

template < class PopulationType >
class CrossoverGenerator {
public:
    typedef PopulationType                                                      population_type;
    typedef typename population_type::allele_space_type::device_space_type      allele_space_type;
//    typedef typename population_type::sequence_space_type                       sequence_space_type;

    typedef typename population_type::real_type             real_type;
    typedef typename population_type::int_type              int_type;
    typedef typename population_type::order_tag_type        order_tag_type;

//    typedef curandState_t                                   state_type;
//    typedef clotho::cuda::curand_helper< state_type >       helper_type;
//    typedef typename helper_type::seed_type                 seed_type;
    typedef clotho::cuda::curand_state_pool                 state_pool_type;

    typedef poisson_cdf< real_type, 32 >                    poisson_type;
    typedef scaled_mean_helper< real_type, order_tag_type > scaled_mean_type;

//    static const int_type   THREADS_PER_BLOCK   = allele_space_type::ALIGNMENT_SIZE;
//    static const int_type   BLOCKS_PER_GRID     = 100;

    CrossoverGenerator( boost::property_tree::ptree & config ) :
//        m_states( NULL )
//        , m_seed( 0 )
        m_recombination_rate( 0.001 )
    {
        parse_configuration( config );

        initialize();
    }

    void operator()( population_type * pop ) {
//        std::cerr << "Generating Random Crossover" << std::endl;
//        crossover_kernel<<< BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>( m_states, pop->alleles.get_device_space(), dPoisCDF, pop->sequences.get_device_space() );
        unsigned int bcount = state_pool_type::getInstance()->get_max_blocks();
        unsigned int tcount = state_pool_type::getInstance()->get_max_threads();
        crossover_kernel<<< bcount, tcount >>>( state_pool_type::getInstance()->get_device_states(), pop->alleles.get_device_space(), pop->free_space, dPoisCDF, pop->sequences.get_device_space() );
    }

    virtual ~CrossoverGenerator() {
//        helper_type::cleanup_states( m_states );
        cudaFree( dPoisCDF );
    }

protected:

    void parse_configuration( boost::property_tree::ptree & config ) {
        boost::property_tree::ptree lconfig;

        if( config.get_child_optional( "recombination" ) != boost::none ) {
            lconfig = config.get_child( "recombination" );
        }

        if( lconfig.get_child_optional( "crossover_per_sequence" ) == boost::none ) {
            lconfig.put("crossover_per_sequence", m_recombination_rate );
        } else {
            m_recombination_rate = lconfig.get< real_type >( "crossover_per_sequence", m_recombination_rate );
        }

//        if( lconfig.get_child_optional( "rng.seed" ) != boost::none ) {
//            m_seed = lconfig.get< seed_type >( "rng.seed", m_seed );
//        }
//
//        if( m_seed == 0 ) {
//            m_seed = clotho::utility::clock_type::now().time_since_epoch().count();
//            lconfig.put("rng.seed", m_seed );
//        }

        state_pool_type::getInstance()->initialize( config );

        config.put_child( "recombination", lconfig );
    }

    void initialize() {

        assert( cudaMalloc( (void **) &dPoisCDF, sizeof( poisson_type) ) == cudaSuccess );

        initialize_poisson( m_recombination_rate, (order_tag_type *) NULL );

//        helper_type::make_states( m_states, m_seed, BLOCKS_PER_GRID, THREADS_PER_BLOCK );
    }

    void initialize_poisson( real_type mean, unordered_tag * p ) {
//        real_type lambda =  ( mean / (real_type) THREADS_PER_BLOCK);
        real_type lambda = (mean / (real_type) allele_space_type::ALIGNMENT_SIZE);

//        std::cerr << "Scaled recombination rate for unordered space: " << mean << " -> " << lambda << std::endl;
        make_poisson_cdf_maxk32<<< 1, 32 >>>( dPoisCDF, lambda );
    }

    void initialize_poisson( real_type mean, unit_ordered_tag< int_type > * p ) {
        real_type lambda = (mean / 32.0);

//        std::cerr << "Scaled recombination rate for unit ordered space: " << mean << " -> " << lambda << std::endl;

        make_poisson_cdf_maxk32<<< 1, 32 >>>( dPoisCDF, lambda );
    }

//    state_type  * m_states;
//    seed_type   m_seed;

    poisson_type * dPoisCDF;
    real_type   m_recombination_rate;
};

#endif  // CROSSOVER_GENERATOR_HPP_
