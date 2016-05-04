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
#include <boost/random/uniform_int_distribution.hpp>

#include "clotho/data_spaces/allele_space/allele_generator.hpp"

namespace clotho {
namespace genetics {

template < class RNG, class GeneticSpaceType >
class MutationGenerator {
public:
    typedef GeneticSpaceType                            genetic_space_type;
    typedef typename genetic_space_type::allele_type    allele_type;

    typedef RNG                                 random_engine_type;

    typedef boost::random::uniform_int_distribution< size_t >   sequence_generator_type;
    typedef AlleleGenerator< random_engine_type, allele_type >  allele_generator_type;

    MutationGenerator( random_engine_type * rng, boost::property_tree::ptree & config ) :
        m_rand( rng )
        , m_allele_gen( rng, config )
    {}

    void operator()( genetic_space_type * child, size_t N ) {
        m_seq_gen.param( typename sequence_generator_type::param_type( 0, child->sequence_count() ) );

        while( N-- ) {
            size_t seq_idx = 1;
            size_t all_idx = child->getAlleleSpace().next_free();

            assert( all_idx != -1 );

            child->getSequenceSpace().flip( seq_idx, all_idx );
            m_allele_gen( child->getAlleleSpace(), all_idx );
        }
    }

    virtual ~MutationGenerator() {}

protected:
    random_engine_type      * m_rand;
    allele_generator_type   m_allele_gen;
    sequence_generator_type m_seq_gen;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_MUTATION_GENERATOR_HPP_
