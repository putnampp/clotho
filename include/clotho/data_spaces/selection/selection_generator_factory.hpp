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
#ifndef CLOTHO_SELECTION_GENERATOR_FACTORY_HPP_
#define CLOTHO_SELECTION_GENERATOR_FACTORY_HPP_

#include "clotho/data_spaces/selection/selection_details.hpp"
#include "clotho/data_spaces/fitness/general_fitness.hpp"

#include "clotho/data_spaces/selection/selection_generator_fitness_inline.hpp"

namespace clotho {
namespace genetics {

template < class RNG, class FitnessSpaceType >
static std::shared_ptr< selection_details< RNG, unsigned int > > make_selection_generator( RNG * rng, FitnessSpaceType * fspace );

template < class RNG >
static std::shared_ptr< selection_details< RNG > > make_selection_generator( RNG * rng, GeneralFitness * fspace  ) {
    return make_selection_generator( rng, fspace, 0, fspace->size() );
}

template < class RNG >
static std::shared_ptr< selection_details< RNG > > make_selection_generator( RNG * rng, GeneralFitness * fspace, unsigned int start, unsigned int end  ) {

    assert( fspace != NULL );


    if( fspace->isEmpty() )  {
        std::shared_ptr< ConstantPairGenerator< RNG > > tmp( new ConstantPairGenerator< RNG >( rng, 0 ) );
        return tmp;
    } else {

        typedef GeneralFitness::iterator iterator;
        bool constant_fitness = true;
        unsigned int N = 0;
        iterator first = fspace->begin() + start;
        iterator last = fspace->begin() + end;

        if( first != last ) {
            typename GeneralFitness::fitness_type prev = *first++;
            while( constant_fitness && ( first != last ) ) {
                constant_fitness = ( prev == *first++);
            }
            N = end - start - 1;
        }

        if( constant_fitness ) {
            std::shared_ptr< UniformPairGenerator< RNG > > tmp( new UniformPairGenerator< RNG >( rng, N ) );

            return tmp;
        }

        std::shared_ptr< DiscretePairGenerator< RNG, typename GeneralFitness::fitness_type > > tmp( new DiscretePairGenerator< RNG, typename GeneralFitness::fitness_type>( rng, fspace->begin() + start, fspace->begin() + end) );
        return tmp;
    }
}

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_SELECTION_GENERATOR_FACTORY_HPP_

