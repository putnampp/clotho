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
#ifndef CLOTHO_GENERAL_FITNESS_HPP_
#define CLOTHO_GENERAL_FITNESS_HPP_

#include <boost/property_tree/ptree.hpp>

#include "clotho/fitness/fitness_toolkit.hpp"

namespace clotho {
namespace genetics {

class GeneralFitness {
public:
    typedef std::shared_ptr< ifitness_generator >           fitness_generator;
    typedef std::shared_ptr< ifitness >                     fitness_operator;

    GeneralFitness( boost::property_tree::ptree & config ) :
        m_fit_gen()
    {

        m_fit_gen = fitness_toolkit::getInstance()->get_tool( config );

        if( !m_fit_gen ) {
            fitness_toolkit::getInstance()->tool_configurations( config );
        }
    }

    template < class PopulationSpaceType, class PhenotypeSpaceType >
    void operator()( PopulationSpaceType * pop, const PhenotypeSpaceType & phenos ) {
        typedef typename PhenotypeSpaceType::phenotype_type         weight_type;
        typedef typename PopulationSpaceType::fitness_score_type    fitness_score;

        size_t N = pop->individual_count();

        if( N == 0 ) return;

        fitness_operator op = m_fit_gen->generate( N );

        weight_type * tmp = phenos.getPhenotypes();
        unsigned int traits = phenos.trait_count();

        unsigned int i = 0;
        while( i < N ) {
            fitness_score score = (*op)(tmp, tmp + traits);

            pop->setFitnessAt( i, score );

            tmp += traits;
            ++i;
        }
    }

    virtual ~GeneralFitness() {}

protected:
    fitness_generator   m_fit_gen;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_GENERAL_FITNESS_HPP_
