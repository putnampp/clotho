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
#ifndef CLOTHO_BATCH_PHENOTYPE_TASK_HPP_
#define CLOTHO_BATCH_PHENOTYPE_TASK_HPP_

#include "clotho/data_spaces/task/task.hpp"

#include "clotho/data_spaces/phenotype_evaluator/phenotype_task.hpp"
#include "clotho/data_spaces/phenotype_evaluator/phenotype_accumulator.hpp"

namespace clotho {
namespace genetics {

template < class SequenceSpaceType, class TraitSpaceType >
class batch_phenotype_task : public task {
public:

    typedef batch_phenotype_task< SequenceSpaceType, TraitSpaceType >   self_type;

    typedef SequenceSpaceType                                           sequence_space_type;
    typedef TraitSpaceType                                              trait_space_type;

    typedef typename sequence_space_type::block_type                    block_type;
    typedef typename trait_space_type::weight_type                      weight_type;
    typedef phenotype_accumulator< trait_space_type >                   accumulator_type;

    typedef phenotype_task< block_type, accumulator_type >              task_type;

    typedef typename sequence_space_type::sequence_vector  sequence_vector;

    batch_phenotype_task( sequence_space_type * pop, trait_space_type * traits, unsigned int start, unsigned int end, weight_type * phenos ) :
        m_pop( pop )
        , m_traits( traits )
        , m_start( start )
        , m_end( end )
        , m_phenos( phenos )
    {    }

    batch_phenotype_task( const self_type & other ) :
        m_pop( other.m_pop )
        , m_traits( other.m_traits )
        , m_start( other.m_start )
        , m_end( other.m_end )
        , m_phenos( other.m_phenos )
    {}

    void operator()() {

        weight_type * tmp = m_phenos;

        for( unsigned int i = m_start; i < m_end; ++i ) {
            accumulator_type acc( m_traits, tmp );

            sequence_vector s = m_pop->getSequence( i );

            task_type t( s.first, s.second, acc );

            t();

            tmp += m_traits->trait_count();
        }
    }

    virtual ~batch_phenotype_task() {}

protected:
    sequence_space_type     * m_pop;
    trait_space_type        * m_traits;

    unsigned int m_start, m_end;

    weight_type             * m_phenos;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_BATCH_PHENOTYPE_TASK_HPP_

