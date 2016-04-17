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
#ifndef CLOTHO_POPCOUNT_ACCUMULATOR_HPP_
#define CLOTHO_POPCOUNT_ACCUMULATOR_HPP_

#include "clotho/utility/popcount.hpp"
#include <vector>

template < class GeneticSpace, class IntType = size_t >
class PopcountAccumulator {
public:

    typedef GeneticSpace genetic_space_type;
    typedef IntType     int_type;

    typedef std::vector< int_type >    accumulator_type;

    PopcountAccumulator( const genetic_space_type & genomes ) {
        update( genomes );
    }

    void update( const genetic_space_type & genomes ) {

        m_counts.reserve( genomes.sequence_count() );

        typedef typename genetic_space_type::const_block_iterator const_block_iterator;

        typedef typename genetic_space_type::block_type block_type;

        const_block_iterator block_iter = genomes.block_iterator();
        size_t i = 0, j = 0;
        while( block_iter.hasNext() ) {
            block_type b = block_iter.next();

            if( j == 0 ) {
                if( j < m_counts.size() ) {
                    m_counts.push_back(0);
                } else {
                    m_counts[i] = 0;
                }
            }

            m_counts[i] = popcount( b );

            if( ++i >= N ) {
                i = 0;
                j += genetic_space_type::bit_helper_type::BITS_PER_BLOCK;
            }
        }
    }

    virtual ~PopcountAccumulator() {}
protected:
    accumulator_type   m_counts;    
};

#endif  // CLOTHO_POPCOUNT_ACCUMULATOR_HPP_
