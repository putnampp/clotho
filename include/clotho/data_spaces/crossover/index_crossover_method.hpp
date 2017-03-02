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
#ifndef CLOTHO_INDEX_CROSSOVER_METHOD_HPP_
#define CLOTHO_INDEX_CROSSOVER_METHOD_HPP_

#include "clotho/data_spaces/population_space/population_spaces.hpp"

namespace clotho {
namespace genetics {

template < class ClassifierType, class PopulationType >
class index_crossover_method;

template < class ClassifierType, class IndexType, class BlockType, class WeightType >
class index_crossover_method< ClassifierType, population_space< index_vector_alignment< IndexType, BlockType >, trait_space_vector< WeightType > > > {
public:
    typedef ClassifierType  classifier_type;
    typedef population_space< index_vector_alignment< IndexType, BlockType >, trait_space_vector< WeightType > > space_type;

    typedef index_crossover_method< classifier_type, space_type > self_type;

    typedef typename space_type::genome_pointer     genome_pointer;
    typedef typename space_type::row_pointer        row_pointer;

    index_crossover_method( const classifier_type & cls ) : m_cfier( cls ) {}

    void operator()( genome_pointer & p0_start, genome_pointer & p0_end, genome_pointer & p1_start, genome_pointer & p1_end, row_pointer & offspring ) {

        while( true ) {
            if( p0_start == p0_end ) {
                while( p1_start != p1_end ) {
                    if( !m_cfier( *p1_start) ) {
                        offspring->push_back(*p1_start);
                    }
                    ++p1_start;
                }
                break;
            } else if( p1_start == p1_end ) {
                while( p0_start != p0_end ) {
                    if( m_cfier( *p0_start) ) {
                        offspring->push_back(*p0_start);
                    }
                    ++p0_start;
                }
                break;
            } else if( *p0_start == *p1_start ) {
                // homozygous
                //
                offspring->push_back(*p0_start);

                ++p0_start;
                ++p1_start;
            } else if( *p0_start < *p1_start ) {
                if( m_cfier( *p0_start ) ) {
                    offspring->push_back(*p0_start);
                }
                ++p0_start;
            } else {
                if( !m_cfier( *p1_start) ) {
                    offspring->push_back(*p1_start);
                }
                ++p1_start;
            }
        }
    }

    virtual ~index_crossover_method() {}

protected:
    classifier_type m_cfier;
};

}   // namespace genetics
}   // namespace clotho
#endif  // CLOTHO_INDEX_CROSSOVER_METHOD_HPP_
