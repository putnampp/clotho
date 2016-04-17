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
#ifndef CLOTHO_MUTATION_GENERATOR_HPP_
#define CLOTHO_MUTATION_GENERATOR_HPP_

#include <boost/property_tree/ptree.hpp>

#include "clotho/data_spaces/allele_space/allele_generator.hpp"

namespace clotho {
namepsace genetics {

template < class GeneticSpaceType >
class MutationGenerator {
public:
    typedef GeneticSpaceType genetic_space_type;
    tyepdef genetic_space_type::allele_type allele_type;

    typedef AlleleGenerator< allele_type > allele_generator_type;

    MutationGenerator( boost::property_tree::ptree & config ) :
        allele_gen( config )
{}

    size_t generateSize( genetic_space_type * parent ) {
        return 1;
    }

    void update( genetic_space_type * child, size_t N ) {
        while( N-- ) {
            size_t seq_idx = 1;
            size_t all_idx = child->getAlleleSpace().next_free();

            assert( all_idx != -1 );

            child->getSequenceSpace().flip( seq_idx, all_idx );
            all_gen.update( child->getAlleleSpace(), all_idx );
        }
    }

    virtual ~MutationGenerator() {}

protected:

    allele_generator_type   allele_gen;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_MUTATION_GENERATOR_HPP_
