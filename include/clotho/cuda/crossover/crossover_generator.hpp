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

#ifndef CROSSOVER_VERSION
#define CROSSOVER_VERSION 1
#endif  // CROSSOVER_VERSION

template < class OrderTag, unsigned char V >
struct kernel_exe_updater {
    void operator()( dim3 & bcount, dim3 & tcount ) { }
};

template < class IntType >
struct kernel_exe_updater< unit_ordered_tag< IntType >, 2 > {
    void operator()( dim3 & bcount, dim3 & tcount ) {
        assert( 2560 <= bcount.x * tcount.x * tcount.y);
        bcount.x = 10;
        tcount.x = 32;
        tcount.y = 8;
    }
};

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

    typedef clotho::utility::algo_version< CROSSOVER_VERSION > algorithm_version_type;

    CrossoverGenerator( boost::property_tree::ptree & config ) :
        m_recombination_rate( config )
        , ver( NULL )
    {
        parse_configuration( config );

        initialize();
    }

    void operator()( population_type * pop ) {
        dim3 bcount(1,1,1), tcount(1,1,1);
        bcount.x = state_pool_type::getInstance()->get_max_blocks();
        tcount.x = state_pool_type::getInstance()->get_max_threads();

        kernel_exe_updater< order_tag_type, CROSSOVER_VERSION > upd;
        upd( bcount, tcount );

        crossover_kernel<<< bcount, tcount >>>( state_pool_type::getInstance()->get_device_states()
                                                , pop->alleles.get_device_space()
                                                , pop->free_space
                                                , dPoisCDF
                                                , pop->sequences.get_device_space()
                                                , ver );
    }

    virtual ~CrossoverGenerator() {
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

    void initialize_poisson( real_type mean, unordered_tag * p ) {
#if CROSSOVER_VERSION == 2
        real_type lambda = (mean / 32.0 );
#else
        real_type lambda = (mean / (real_type) allele_space_type::ALIGNMENT_SIZE);
#endif  // CROSSOVER_VERSION

        make_poisson_cdf_maxk32<<< 1, 32 >>>( dPoisCDF, lambda );
    }

    void initialize_poisson( real_type mean, unit_ordered_tag< int_type > * p ) {
        real_type lambda = (mean / 32.0);

        make_poisson_cdf_maxk32<<< 1, 32 >>>( dPoisCDF, lambda );
    }

    poisson_type * dPoisCDF;
    recombination_rate_type  m_recombination_rate;
    algorithm_version_type      * ver;
};

#endif  // CROSSOVER_GENERATOR_HPP_
