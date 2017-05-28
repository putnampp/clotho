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
#ifndef FIT_SELECTION_GENERATOR_HOST_HPP_
#define FIT_SELECTION_GENERATOR_HOST_HPP_

#include "clotho/cuda/curand_state_pool.hpp"

#include "clotho/cuda/data_spaces/basic_data_space.hpp"

#include "clotho/cuda/recombination/recombine_parents.hpp"
#include "clotho/recombination/sequence_bias_parameter.hpp"

#include "clotho/cuda/device_state_object.hpp"

#include "clotho/random/seed_parameter.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/bernoulli_distribution.hpp>

#include "clotho/cuda/selection/fit_selection_generator_def.hpp"

#include "clotho/utility/bit_helper.hpp"

struct host_discrete_distribution {};

template < class IntType, class RealType >
class FitSelectionGenerator< IntType, RealType, host_discrete_distribution > : public clotho::utility::iStateObject {
public:
//    typedef basic_data_space< IntType >                     event_space_type;
    typedef IntType                                       event_space_type;

//    typedef clotho::cuda::curand_state_pool             state_pool_type;
//
//    typedef DiscreteDistribution< IntType, RealType >   discrete_dist_type;
    typedef boost::random::mt19937                                       random_engine_type;
    typedef boost::random::discrete_distribution< IntType, RealType >   discrete_dist_type;
    typedef boost::random::bernoulli_distribution< double >             bernoulli_type;
    typedef sequence_bias_parameter< RealType >         sequence_bias_type;

    typedef clotho::utility::BitHelper< IntType >      bit_helper_type;

    FitSelectionGenerator( boost::property_tree::ptree & config ) :
         m_seq_bias( config )
        , dParentIndices(NULL)
        , localParentIndices( NULL )
        , localFitness( NULL )
        , m_parent_size( 0 )
        , m_current_parent_size(0)
        , m_fitness_capacity(0)
    {
        seed_parameter< > seed( config );

        m_rng.seed( seed.m_seed );

        m_bern.param( m_seq_bias.m_bias );
    }

    void generate( RealType * fitness, unsigned int N ) {

        buildParentIndices( fitness, N );

        cudaError_t err = cudaMemcpy(dParentIndices, localParentIndices, N * sizeof( IntType ), cudaMemcpyHostToDevice );
        if( err != cudaSuccess ) {
            std::cerr << "ERROR: " << cudaGetErrorString( err ) << std::endl;
            assert( false );
        }

    }

    template < class PopulationSpaceType >
    void generate_and_recombine( PopulationSpaceType * parent_pop, PopulationSpaceType * child_pop, RealType * fitness, unsigned int N ) {
        generate( fitness, N );

        clotho::utility::algo_version< 3 > * v = NULL;
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
        boost::property_tree::ptree evts;

        for( unsigned int i = 0; i < m_parent_size; ++i ) {
            clotho::utility::add_value_array( evts, localParentIndices[i] );
        }

        state.add_child( "events", evts );
    }

    virtual ~FitSelectionGenerator() {
        if( localParentIndices != NULL ) {
            delete [] localParentIndices;
            cudaFree( dParentIndices );
        }

        if( localFitness != NULL ) {
            delete [] localFitness;
        }
    }

protected:

    void resize( unsigned int N ) {
        if( N > m_parent_size ) {
            if( localParentIndices != NULL ) {
                delete [] localParentIndices;
                cudaFree( dParentIndices );
            }

            assert( cudaMalloc( (void **) &dParentIndices, N * sizeof(IntType ) ) == cudaSuccess );

            localParentIndices = new IntType[ N ];
            m_parent_size = N;
        }

        if( N > m_fitness_capacity ) {
            if( localFitness != NULL ) {
                delete [] localFitness;
            }

            localFitness = new RealType[ N ];

            m_fitness_capacity = N;
        }
    }

    void buildParentIndices( RealType * fitness, unsigned int N ) {
        resize( N );

        if( m_current_parent_size == 0 ) {
            for( unsigned int i = 0; i < N; ++i ) {
                localParentIndices[ i ] = 0;
            }
        } else {

            cudaError_t err = cudaMemcpy( localFitness, fitness, m_current_parent_size * sizeof( RealType ), cudaMemcpyDeviceToHost );
            if( err != cudaSuccess ) {
                std::cerr << "ERROR: " << cudaGetErrorString( err ) << std::endl;
                assert(false);
            }

            discrete_dist_type dist( localFitness, localFitness + m_current_parent_size );

            for( unsigned int i = 0; i < N; ++i ) {
               
                // sequence offset
                localParentIndices[ i ] = 2 * dist( m_rng );

                if( m_bern( m_rng ) ) {
                    localParentIndices[ i ] += 1;  
                }
            }
        }

        m_current_parent_size = N / 2;
    }

    random_engine_type  m_rng;
    sequence_bias_type  m_seq_bias;
    bernoulli_type  m_bern;

    IntType         * dParentIndices;
    IntType         * localParentIndices;

    RealType        * localFitness;

    unsigned int m_parent_size;
    unsigned int m_current_parent_size, m_fitness_capacity;
};

#endif  // FIT_SELECTION_GENERATOR_HOST_HPP_
