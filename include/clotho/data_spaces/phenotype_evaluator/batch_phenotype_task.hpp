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

template < class GeneticSpaceType >
class batch_phenotype_task : public task {
public:

    typedef batch_phenotype_task< GeneticSpaceType >    self_type;

    typedef GeneticSpaceType                                genetic_space_type;

    typedef typename genetic_space_type::block_type                 block_type;
    typedef typename genetic_space_type::allele_type::weight_type   weight_type;
    typedef phenotype_accumulator< weight_type >                    accumulator_type;

    typedef phenotype_task< block_type, accumulator_type >                  task_type;

    typedef typename genetic_space_type::association_type::sequence_vector  sequence_vector;

    batch_phenotype_task( genetic_space_type * pop, unsigned int start, unsigned int end, unsigned int traits, weight_type * weights  ) :
        m_pop( pop )
        , m_start( start )
        , m_end( end )
        , m_trait_count( traits )
        , m_weights( weights )
    {    }

    batch_phenotype_task( const self_type & other ) :
        m_pop( other.m_pop )
        , m_start( other.m_start )
        , m_end( other.m_end )
        , m_trait_count( other.m_trait_count )
        , m_weights( other.m_weights )
    {}

    void operator()() {

        weight_type * tmp = m_weights;

        for( unsigned int i = m_start; i < m_end; ++i ) {
            accumulator_type acc( m_pop->getAlleleSpace().getWeights(), m_pop->getAlleleSpace().allele_count(), m_trait_count, tmp );

            sequence_vector s = m_pop->getSequenceSpace().getSequence( i );

            task_type t( s.first, s.second, acc );

            t();

            tmp += m_trait_count;
        }
    }

    virtual ~batch_phenotype_task() {}

protected:
    genetic_space_type  * m_pop;
    unsigned int m_start, m_end, m_trait_count;

    weight_type         * m_weights;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_BATCH_PHENOTYPE_TASK_HPP_

