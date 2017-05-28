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
#ifndef HOST_FITNESS_TRANSLATOR_HPP_
#define HOST_FITNESS_TRANSLATOR_HPP_

#include "clotho/mutation/mutation_rate_parameter.hpp"
#include "clotho/fitness/quadratic_fitness_parameter.hpp"

#include "clotho/cuda/fitness/fitness_evaluators.hpp"

template < class RealType >
class HostFitnessTranslator {
public:

    HostFitnessTranslator( boost::property_tree::ptree & config ) :
        m_quad_param( config )
        , m_mutate_rate( config )
    {}

    template < class IntType >
    void operator()( HostPopulationSpace< RealType, IntType > * pop ) {

        dim3 blocks( pop->getIndividualCount(), 1, 1), threads( 1,1,1 );

        RealType stdev = (RealType) pop->getSequenceCount();
        stdev *= 2.0 * m_mutate_rate.m_mu;
        stdev *= m_quad_param.m_scale;

        evaluate_quadritic_fitness<<< blocks, threads >>>( pop->getDevicePhenotypes(), pop->getDeviceFitness(), pop->getTraitCount(), stdev );
    }

    virtual ~HostFitnessTranslator() {}

protected:

    quadratic_fitness_parameter< RealType > m_quad_param;
    mutation_rate_parameter< RealType >     m_mutate_rate;
};

#endif  // HOST_FITNESS_TRANSLATOR_HPP_
