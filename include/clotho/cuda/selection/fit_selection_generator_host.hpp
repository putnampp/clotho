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
    typedef basic_data_space< IntType >                     event_space_type;

//    typedef clotho::cuda::curand_state_pool             state_pool_type;
//
//    typedef DiscreteDistribution< IntType, RealType >   discrete_dist_type;
    typedef boost::random::mt19937                                       random_engine_type;
    typedef boost::random::discrete_distribution< IntType, RealType >   discrete_dist_type;
    typedef boost::random::bernoulli_distribution< double >             bernoulli_type;
    typedef sequence_bias_parameter< RealType >         sequence_bias_type;

    typedef clotho::utility::BitHelper< IntType >      bit_helper_type;

    FitSelectionGenerator( boost::property_tree::ptree & config ) :
        dParentIndices(NULL)
        , m_seq_bias( config )
        , localParentIndices( NULL )
        , m_parent_size( 0 )
    {
//        state_pool_type::getInstance()->initialize( config );
        create_space( dParentIndices );

        seed_parameter< > seed( config );

        m_rng.seed( seed.m_seed );
    }

    template < class PopulationSpaceType, class FitnessSpaceType >
    void generate( PopulationSpaceType * parent_pop, PopulationSpaceType * child_pop, FitnessSpaceType * fitness, unsigned int N ) {
        FitnessSpaceType local_fit;
        cudaError_t err = cudaMemcpy( &local_fit, fitness, sizeof( FitnessSpaceType ), cudaMemcpyDeviceToHost );
        if( err != cudaSuccess ) {
            std::cerr << "ERROR: " << cudaGetErrorString( err ) << std::endl;
            assert(false);
        }

        buildDiscrete( local_fit, N );

        // resize the devices vector
        resize_space( dParentIndices, N );
        // update the hosts copy of the device object (gain device pointer)
        err = cudaMemcpy( &hParentIndices, dParentIndices, sizeof( event_space_type ), cudaMemcpyDeviceToHost );
        if( err != cudaSuccess ) {
            std::cerr << "ERROR: " << cudaGetErrorString( err ) << std::endl;
            assert( false );
        }

        err = cudaMemcpy( hParentIndices.data, localParentIndices, N * sizeof( IntType ), cudaMemcpyHostToDevice );
        if( err != cudaSuccess ) {
            std::cerr << "ERROR: " << cudaGetErrorString( err ) << std::endl;
            assert( false );
        }
    }

    template < class PopulationSpaceType, class FitnessSpaceType >
    void generate_and_recombine( PopulationSpaceType * parent_pop, PopulationSpaceType * child_pop, FitnessSpaceType * fitness, unsigned int N ) {
        generate( parent_pop, child_pop, fitness, N );

        clotho::utility::algo_version< 2 > * v = NULL;
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
//        boost::property_tree::ptree fit;
//        m_discrete.get_state( fit );

        boost::property_tree::ptree evts;
        get_device_object_state( evts, dParentIndices );

//        state.add_child( "distribution", fit );
        state.add_child( "events", evts );
    }

    virtual ~FitSelectionGenerator() {
        delete_space( dParentIndices );
    }

protected:

    void resize( unsigned int N ) {
        if( N > m_parent_size ) {
            if( localParentIndices != NULL ) {
                delete [] localParentIndices;
            }

            localParentIndices = new IntType[ N ];

            m_parent_size = N;
        }
    }

    template < class WeightType >
    void buildDiscrete( basic_data_space< WeightType > & discrete_weights, unsigned int N ) {
        resize( N );

        typedef typename basic_data_space< WeightType >::value_type value_type;

        value_type * data = new value_type[ discrete_weights.size ];

        copy_heap_data( data, discrete_weights.data, discrete_weights.size );

        discrete_dist_type dist( data, data + discrete_weights.size );

        for( unsigned int i = 0; i < N; ++i ) {
           
            localParentIndices[ i ] = dist( m_rng );

            if( m_bern( m_rng ) ) {
                localParentIndices[ i ] |= bit_helper_type::MSB_SET;  
            }
        }

        delete [] data;
    }

    random_engine_type  m_rng;
    event_space_type    hParentIndices;
    event_space_type    * dParentIndices;

    sequence_bias_type  m_seq_bias;

    bernoulli_type  m_bern;
    IntType         * localParentIndices;

    RealType        * localWeights;

    unsigned int m_parent_size;
};

#endif  // FIT_SELECTION_GENERATOR_HOST_HPP_
