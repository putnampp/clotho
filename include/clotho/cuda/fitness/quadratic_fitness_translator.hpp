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
#ifndef QUADRATIC_FITNESS_TRANSLATOR_HPP_
#define QUADRATIC_FITNESS_TRANSLATOR_HPP_

#include "clotho/fitness/fitness_translator.hpp"

#include "clotho/cuda/data_spaces/basic_data_space.hpp"
#include "clotho/mutation/mutation_rate_parameter.hpp"

#include "clotho/cuda/device_state_object.hpp"

#include <cuda.h>

template < class RealType >
__global__ void quadratic_fitness_kernel( basic_data_space< RealType > * phenos, RealType scale_coeff, basic_data_space< RealType > * fitness ) {
    typedef RealType real_type;

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    unsigned int N = phenos->size;
    unsigned int M = fitness->size;

    assert( N <= M );

    real_type * pdata = phenos->data;
    real_type * fdata = fitness->data;

    while( tid < N ) {
        real_type x = pdata[tid];

        x /= scale_coeff;
        x *= x;

        // x = 1.0 - ((real_type)(x > 1.0)* 1.0) - ((real_type)(x < 1.0) * x)
        x = (( x > 1.0 ) ? 0.0 : (1.0 - x));

        fdata[tid] = x;
        tid += (blockDim.x * blockDim.y);
    }
}

template < class PhenotypeSpaceType >
class QuadraticFitnessTranslator :
    public fitness_translator< PhenotypeSpaceType >
{
public:
    typedef PhenotypeSpaceType                          phenotype_space_type;
    typedef typename phenotype_space_type::value_type   real_type;

    typedef basic_data_space< real_type >        fitness_space_type;

    QuadraticFitnessTranslator( boost::property_tree::ptree & config ) :
        dFitness( NULL )
        , m_mutation_rate( config )
        , m_scale( 1.0 )
    {
        parse_configuration( config );
        create_space( dFitness );
    }

    void operator()( phenotype_space_type * phenos, unsigned int N ) {
        resize_space( dFitness, N / 2 );

        real_type stdev = 2.0 * ((real_type)N) * m_mutation_rate.m_mu;
        stdev = sqrt( stdev );

        stdev *= m_scale;

        quadratic_fitness_kernel<<< 1, 1024 >>>( phenos, stdev, dFitness );
    }

    fitness_space_type * get_device_space() {
        return dFitness;
    }

    void get_state( boost::property_tree::ptree & state ) {
        get_device_object_state( state, dFitness );
    }

    virtual ~QuadraticFitnessTranslator() {
        delete_space( dFitness );
    }

protected:

    void parse_configuration( boost::property_tree::ptree & config ) {
        boost::property_tree::ptree lconfig;

        if( config.get_child_optional( "fitness" ) != boost::none ) {
            lconfig = config.get_child( "fitness" );
        }

        if( lconfig.get_child_optional( "scale" ) == boost::none ) {
            lconfig.put( "scale", m_scale );
        } else {
            m_scale = lconfig.get< real_type >( "scale", m_scale );

            assert( m_scale > 0.0 );
        }

        config.put_child( "fitness", lconfig );
    }

    fitness_space_type  * dFitness;

    mutation_rate_parameter< real_type > m_mutation_rate;
    real_type   m_scale;
};

#endif  // QUADRATIC_FITNESS_TRANSLATOR_HPP_
