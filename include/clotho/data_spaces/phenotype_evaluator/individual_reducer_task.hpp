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
#ifndef CLOTHO_INDIVIDUAL_REDUCER_HPP_
#define CLOTHO_INDIVIDUAL_REDUCER_HPP_

#include "clotho/data_spaces/task/task.hpp"

namespace clotho {
namespace genetics {

template < class SpaceType, class PhenotypeSpaceType >
class individual_reducer_task : public task {
public:
    typedef individual_reducer_task< SpaceType, PhenotypeSpaceType >                    self_type;
    typedef SpaceType                                                                   space_type;

    typedef PhenotypeSpaceType                                                          phenotype_space_type;

    typedef typename space_type::base_genome_type::weight_type::weight_vector           weight_vector;
    typedef typename space_type::base_genome_type::weight_type::const_weight_iterator   weight_iterator;

    typedef typename space_type::individual_type                                        individual_type;
    typedef typename space_type::individual_iterator                                    individual_iterator;
    typedef typename space_type::const_individual_iterator                              const_individual_iterator;

    typedef typename space_type::genome_type                                            genome_type;

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
        unsigned int i = 0;

        individual_iterator first = m_pop->begin_individual();

        while( i < m_start ) { ++first; ++i; }

        while( i < m_end ) {
            weight_vector w = m_pop->create_weight_vector();
            
            combine( first->first, first->second, w );

            ++first;
            ++i;
        }
    }

    virtual ~individual_reducer_task() {}

protected:

    void combine( const genome_type & s0, const genome_type & s1, weight_vector & res ) {
        weight_iterator s0_first = s0->begin_weight(), s0_last = s0->end_weight();
        weight_iterator s1_first = s1->begin_weight(), s1_last = s1->end_weight();

        while( true ) {
            if( s0_first == s0_last ) {
#ifdef DEBUGGING
                assert( s1_first == s1_last );
#endif  // DEBUGGING
                break;
            } else if( s1_first == s1_last ) {
#ifdef DEBUGGING
                assert( false );
#endif  // DEBUGGING
                break;
            }

            res.push_back( *s0_first + *s1_first );

            ++s0_first;
            ++s1_first;
        }
    }

    space_type  * m_pop;

    unsigned int    m_start, m_end;

    phenotype_space_type * m_phenos;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_INDIVIDUAL_REDUCER_HPP_
