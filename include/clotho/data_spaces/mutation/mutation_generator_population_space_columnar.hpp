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
#ifndef CLOTHO_MUTATION_GENERATOR2_POPULATION_SPACE_COLUMNAR_HPP_
#define CLOTHO_MUTATION_GENERATOR2_POPULATION_SPACE_COLUMNAR_HPP_

#include "clotho/data_spaces/population_space/population_space_columnar.hpp"

#include <boost/property_tree/ptree.hpp>
#include <boost/random/uniform_int_distribution.hpp>

namespace clotho {
namespace genetics {

template < class RNG, class BlockType, class WeightType >
class MutationGenerator2 < RNG, population_space_columnar< BlockType, WeightType > > {
public:
    typedef RNG     random_number_engine;

    typedef population_space_columnar< BlockType, WeightType >   sequence_space_type;

    typedef boost::random::uniform_int_distribution< unsigned int > sequence_distribution_type;

    MutationGenerator2( random_number_engine * rng, boost::property_tree::ptree & config ) :
        m_rng(rng)
    {}

    template < class FreeSpaceType >
    void operator()( sequence_space_type * seqs, FreeSpaceType & free_space, unsigned int N  ) {
        assert( N <= free_space.free_size() );

        sequence_distribution_type seq_gen( 0, seqs->individual_count() - 1 );

        typename FreeSpaceType::base_type::iterator it = free_space.free_begin();

        while( N-- ) {
            unsigned int seq_idx = seq_gen( *m_rng );
            typename FreeSpaceType::size_type all_idx = *it++;

            operator()( seqs, seq_idx, all_idx );
        }
    }

    inline void operator()( sequence_space_type * seqs, unsigned int seq_idx, unsigned int all_idx ) {
        //std::cerr << "Mutating " << all_idx << std::endl;
        seqs->mutate( seq_idx, all_idx );
    }

    virtual ~MutationGenerator2() {}

protected:
    random_number_engine * m_rng;
};

}   // namespace genetics
}   // namespace clotho
#endif  // CLOTHO_MUTATION_GENERATOR2_POPULATION_SPACE_COLUMNAR_HPP_
