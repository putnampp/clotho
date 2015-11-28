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
#include "clotho/cuda/helpers/scaled_mean_helper.hpp"

#include "clotho/cuda/curand_state_pool.hpp"

#include "clotho/cuda/distributions/poisson_distribution.hpp"
#include "clotho/recombination/recombination_rate_parameter.hpp"

#ifndef CROSSOVER_VERSION
#define CROSSOVER_VERSION 1
#endif  // CROSSOVER_VERSION

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
    typedef kernel_exec< xover_config< order_tag_type, CROSSOVER_VERSION > > exec_type;

    CrossoverGenerator( boost::property_tree::ptree & config ) :
        m_recombination_rate( config )
        , ver( NULL )
        , m_exec( config )
        , m_blocks( 1, 1, 1)
        , m_threads( 1, 1, 1 )
    {
        parse_configuration( config );

        initialize();
    }

    void operator()( population_type * pop ) {
        CHECK_LAST_KERNEL_EXEC

//        std::cout << "Blocks: {" << m_blocks.x << ", " << m_blocks.y << ", " << m_blocks.z << "}" << std::endl;
//        std::cout << "Threads: {" << m_threads.x << ", " << m_threads.y << ", " << m_threads.z << "}" << std::endl;

        crossover_kernel<<< m_blocks, m_threads >>>( state_pool_type::getInstance()->get_device_states()
                                                , pop->alleles.get_device_space()
                                                , pop->free_space
                                                , dPoisCDF
                                                , pop->sequences.get_device_space()
                                                , ver );
        CHECK_LAST_KERNEL_EXEC
    }

    virtual ~CrossoverGenerator() {
        cudaFree( dPoisCDF );
    }

protected:

    void parse_configuration( boost::property_tree::ptree & config ) {
        state_pool_type::getInstance()->initialize( config );

        m_blocks.x = state_pool_type::getInstance()->get_max_blocks();
        m_threads.x = state_pool_type::getInstance()->get_max_threads();

        m_exec( m_blocks, m_threads );
    }

    void initialize() {

        assert( cudaMalloc( (void **) &dPoisCDF, sizeof( poisson_type) ) == cudaSuccess );


        real_type lambda = compute_lambda( m_recombination_rate.m_rho, (order_tag_type *) NULL, (clotho::utility::algo_version< CROSSOVER_VERSION > *) NULL );

        
        make_poisson_cdf_maxk32<<< 1, 32 >>>( dPoisCDF, lambda );
    }

    template < class Tag0, class Tag1 >
    real_type compute_lambda( real_type mean, Tag0 * t0, Tag1 * t1 ) {
        return mean;
    }

    template < unsigned char V >
    real_type compute_lambda( real_type mean, unordered_tag * p, clotho::utility::algo_version< V > * v ) {
        std::cerr << "Scaling mean by 1/" << allele_space_type::ALIGNMENT_SIZE << std::endl;
        return (mean / (real_type) allele_space_type::ALIGNMENT_SIZE);
    }

    real_type compute_lambda( real_type mean, unordered_tag * p, clotho::utility::algo_version< 2 > * v ) {
        return (mean / 32.0);
    }

    real_type compute_lambda( real_type mean, unordered_tag * p, clotho::utility::algo_version< 3 > * v ) {
        return (mean / 32.0);
    }

    template < unsigned char V >
    real_type compute_lambda( real_type mean, unit_ordered_tag< int_type > * p, clotho::utility::algo_version< V > * ver ) {
        return (mean / 32.0);
    }

    poisson_type * dPoisCDF;
    recombination_rate_type  m_recombination_rate;
    algorithm_version_type      * ver;
    exec_type                   m_exec;
    dim3                        m_blocks, m_threads;
};

#endif  // CROSSOVER_GENERATOR_HPP_
