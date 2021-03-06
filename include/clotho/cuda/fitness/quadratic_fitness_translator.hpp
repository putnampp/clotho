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
#include "clotho/fitness/quadratic_fitness_parameter.hpp"
#include "clotho/mutation/mutation_rate_parameter.hpp"

#include "clotho/cuda/device_state_object.hpp"
#include "clotho/cuda/fitness/quadratic_fitness_kernel.hpp"

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
        , m_quad( config )
        , m_mutation_rate( config )
    {
        create_space( dFitness );
    }

    void operator()( phenotype_space_type * phenos, unsigned int N ) {
        resize_space( dFitness, N / 2 );

        real_type stdev = 2.0 * ((real_type)N) * m_mutation_rate.m_mu;
        stdev = sqrt( stdev );

        stdev *= m_quad.m_scale;

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
    fitness_space_type  * dFitness;

    quadratic_fitness_parameter< real_type > m_quad;
    mutation_rate_parameter< real_type > m_mutation_rate;
};

#endif  // QUADRATIC_FITNESS_TRANSLATOR_HPP_
