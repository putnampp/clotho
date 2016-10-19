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
#ifndef CLOTHO_INDIVIDUAL_REDUCER_COLUMNAR_HPP_
#define CLOTHO_INDIVIDUAL_REDUCER_COLUMNAR_HPP_

#include "clotho/data_spaces/task/task.hpp"

namespace clotho {
namespace genetics {

template < class BlockType, class WeightType, class PhenotypeSpaceType >
class individual_reducer_task< population_space_columnar< BlockType, WeightType >, PhenotypeSpaceType > : public task {
public:
    typedef individual_reducer_task< population_space_columnar< BlockType, WeightType >, PhenotypeSpaceType >                    self_type;
    typedef population_space_columnar< BlockType, WeightType >                          space_type;

    typedef PhenotypeSpaceType                                                          phenotype_space_type;

    typedef typename space_type::weight_vector                                          weight_vector;
    typedef typename space_type::weight_iterator                                        weight_iterator;
    typedef typename space_type::const_weight_iterator                                  const_weight_iterator;

    typedef typename phenotype_space_type::phenotype_iterator                           phenotype_iterator;


    individual_reducer_task( space_type * pop, unsigned int start, unsigned int end, phenotype_space_type * phenos  ) :
        m_pop( pop )
        , m_start(start)
        , m_end(end)
        , m_phenos( phenos )
    {}

    individual_reducer_task( const self_type & other ) :
        m_pop( other.m_pop )
        , m_start( other.m_start )
        , m_end( other.m_end )
        , m_phenos( other.m_phenos )
    {}

    void operator()() {

        for( unsigned int i = m_start; i < m_end; ++i ) {
            phenotype_iterator res = m_phenos->begin_individual_phenotype( i );

            unsigned int j = 2 * i;
            const_weight_iterator s0f = m_pop->begin_genome_traits( j ), s0l = m_pop->end_genome_traits( j );
            ++j;

            const_weight_iterator s1f = m_pop->begin_genome_traits( j ), s1l = m_pop->end_genome_traits( j );

            combine( res, s0f, s0l, s1f, s1l );
        }
    }

    virtual ~individual_reducer_task() {}

protected:

    void combine( phenotype_iterator & res, const_weight_iterator & s0_first, const_weight_iterator & s0_last, const_weight_iterator & s1_first, const_weight_iterator & s1_last ) {

        while( s0_first != s0_last ) {
            *res++ = (*s0_first + *s1_first);
            ++s0_first;
            ++s1_first;
        }

        assert( s1_first == s1_last );
    }


    space_type  * m_pop;

    unsigned int    m_start, m_end;

    phenotype_space_type * m_phenos;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_INDIVIDUAL_REDUCER_COLUMNAR_HPP_
