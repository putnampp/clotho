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
#ifndef CLOTHO_PHENOTYPE_EVALUATOR_FITNESS_HPP_
#define CLOTHO_PHENOTYPE_EVALUATOR_FITNESS_HPP_

#include <boost/property_tree/ptree.hpp>
#include "clotho/data_spaces/phenotype_evaluator/phenotype_evaluator.hpp"

#include "clotho/fitness/fitness_toolkit.hpp"

namespace clotho {
namespace genetics {

template < class TraitAccumType >
class Fitness < phenotype_evaluator< TraitAccumType > > {
public:
    typedef phenotype_evaluator< TraitAccumType >           evaluator_type;
    typedef TraitAccumType                                  trait_accumulator_type;

    typedef typename evaluator_type::phenotype_type         phenotype_type;

    typedef typename evaluator_type::genetic_space_type     genetic_space_type;
    typedef double                                          result_type;

    typedef std::shared_ptr< ifitness_generator >           fitness_generator;
    typedef std::shared_ptr< ifitness >                     fitness_operator;

    Fitness( boost::property_tree::ptree & config ) :
        m_fitness()
    {
        m_fitness = fitness_toolkit::getInstance()->get_tool( config );

        if( !m_fitness ) {
            fitness_toolkit::getInstance()->tool_configurations( config );
        }
    }

    void update( genetic_space_type * pop, evaluator_type & eval ) {
        fitness_operator op = m_fitness->generate( eval.getPhenotypes() );
    }

    virtual ~Fitness() {}

protected:
    fitness_generator    m_fitness;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_PHENOTYPE_EVALUATOR_FITNESS_HPP_

