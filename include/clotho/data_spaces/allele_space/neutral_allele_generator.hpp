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
#ifndef CLOTHO_NEUTRAL_ALLELE_GENERATOR_HPP_
#define CLOTHO_NEUTRAL_ALLELE_GENERATOR_HPP_

#include "clotho/data_spaces/allele_space/neutral_allele.hpp"
#include "clotho/data_spaces/allele_space/allele_generator.hpp"

#include "clotho/data_spaces/generators/position_generator.hpp"
#include "clotho/data_spaces/generators/neutral_generator.hpp"

namespace clotho {
namespace genetics {

template < class RNG, class PositionType >
class AlleleGenerator< RNG, neutral_allele_vectorized< PositionType > > {
public:
    typedef neutral_allele_vectorized< PositionType >   allele_type;

    typedef position_generator< RNG, PositionType >     position_generator_type;
    typedef neutral_generator< RNG >                    neutral_generator_type;

    AlleleGenerator( RNG * rng, boost::property_tree::ptree & config ):
        m_pos_gen( rng, config )
        , m_neut_gen( rng, config )
    {}

    void operator()( allele_type & all, size_t offset ) {
        generate( all, offset );
    }

    void generate( allele_type & all, size_t offset ) {
        all.setPositionAt( offset, m_pos_gen() );
        all.setNeutralAt( offset, m_neut_gen() );
    }

    virtual ~AlleleGenerator() {}
protected:

    position_generator_type     m_pos_gen;
    neutral_generator_type      m_neut_gen;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_NEUTRAL_ALLELE_GENERATOR_HPP_
