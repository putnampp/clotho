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

namespace clotho {
namespace genetics {

template < class TraitAccumType >
class Fitness < phenotype_evaluator< TraitAccumType > > {
public:
    typedef phenotype_evaluator< TraitAccumType >           phenotype_type;
    typedef TraitAccumType                                  trait_accumulator_type;

    typedef typename phenotype_type::genetic_space_type     genetic_space_type;
    typedef double                                          result_type;

    Fitness( boost::property_tree::ptree & config ) {}

    void update( genetic_space_type * pop, phenotype_type & phenos ) {

    }

    result_type getFitnessAt( size_t idx ) {
        return 0.0;
    }

    virtual ~Fitness() {}
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_PHENOTYPE_EVALUATOR_FITNESS_HPP_

