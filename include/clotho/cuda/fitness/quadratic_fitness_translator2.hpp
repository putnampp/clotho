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
#ifndef QUADRATIC_FITNESS_TRANSLATOR_2_HPP_
#define QUADRATIC_FITNESS_TRANSLATOR_2_HPP_

#include "clotho/fitness/fitness_translator.hpp"

#include "clotho/cuda/data_spaces/basic_data_space.hpp"
#include "clotho/fitness/quadratic_fitness_parameter.hpp"
#include "clotho/mutation/mutation_rate_parameter.hpp"

#include "clotho/cuda/device_state_object.hpp"
#include "clotho/cuda/fitness/quadratic_fitness_kernel.hpp"

template < class PhenotypeSpaceType >
class QuadraticFitnessTranslator2 :
    public fitness_translator< PhenotypeSpaceType >
{
public:
    typedef PhenotypeSpaceType                          phenotype_space_type;
    typedef typename phenotype_space_type::value_type   real_type;

//    typedef basic_data_space< real_type >        fitness_space_type;
    typedef real_type                              fitness_space_type;

    QuadraticFitnessTranslator2( boost::property_tree::ptree & config ) :
        dFitness( NULL )
        , m_quad( config )
        , m_mutation_rate( config )
        , m_size(0)
        , m_capacity( 0 )
    { }

    void operator()( phenotype_space_type * phenos, unsigned int N ) {
        resize( N / 2 );

        real_type stdev = 2.0 * ((real_type)N) * m_mutation_rate.m_mu;
        stdev = sqrt( stdev );

        stdev *= m_quad.m_scale;

        quadratic_fitness_kernel<<< 1, 1024 >>>( phenos, stdev, dFitness, N / 2 );

    }

    fitness_space_type * get_device_space() {
        return dFitness;
    }

    void get_state( boost::property_tree::ptree & state ) {
    //    get_device_object_state( state, dFitness );
    //
        fitness_space_type * local = new fitness_space_type[ m_capacity ];

        assert( cudaMemcpy( local, dFitness, m_capacity * sizeof( fitness_space_type ), cudaMemcpyDeviceToHost ) == cudaSuccess );

        boost::property_tree::ptree tmp;
        for( unsigned int i = 0; i < m_size; ++i ) {
            clotho::utility::add_value_array(tmp, local[ i ] );
        }

        state.put( "size", m_size );
        state.put( "capacity", m_capacity );
        state.add_child("data", tmp );

        delete [] local;
    }

    virtual ~QuadraticFitnessTranslator2() {

        if( dFitness != NULL ) {
            cudaFree( dFitness );
        }
    }

protected:

    void resize( unsigned int N ) {
        if( N > m_capacity ) {
            if( dFitness != NULL ) {
                cudaFree( dFitness );
            }

            assert( cudaMalloc( (void **) &dFitness, N * sizeof( fitness_space_type ) ) == cudaSuccess );

            m_capacity = N;
        }

        m_size = N;
    }

    fitness_space_type  * dFitness;

    quadratic_fitness_parameter< real_type > m_quad;
    mutation_rate_parameter< real_type > m_mutation_rate;

    unsigned int m_size, m_capacity;
};

#endif  // QUADRATIC_FITNESS_TRANSLATOR_2_HPP_
