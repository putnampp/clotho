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

    typedef typename genetic_space_type::fitness_score_type fitness_score;

    Fitness( boost::property_tree::ptree & config ) :
        m_fit_gen()
    {
        m_fit_gen = fitness_toolkit::getInstance()->get_tool( config );

        if( !m_fit_gen ) {
            fitness_toolkit::getInstance()->tool_configurations( config );
        }
    }

    void update( genetic_space_type * pop, evaluator_type & eval ) {
        size_t N = eval.getPhenotypes().size();

        if( N == 0 )    return;

//        std::cerr << N << " == " << pop->individual_count() << " [ " << pop->sequence_count() << " ]"  << std::endl;
        assert( N == pop->individual_count() );
        fitness_operator op = m_fit_gen->generate( N );

        size_t i = 0;

        while( i < N ) {
            fitness_score f = (*op)( eval.getPhenotypeAt(i) );
            pop->setFitnessAt( i, f);
            ++i;
        }
    }

    virtual ~Fitness() {}

protected:
    fitness_generator   m_fit_gen;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_PHENOTYPE_EVALUATOR_FITNESS_HPP_

