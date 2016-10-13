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
#ifndef CLOTHO_BATCH_PHENOTYPE_TASK_POPULATION_SPACE_HPP_
#define CLOTHO_BATCH_PHENOTYPE_TASK_POPULATION_SPACE_HPP_

#include "clotho/data_spaces/task/task.hpp"

#include "clotho/utility/debruijn_bit_walker.hpp"

#include "clotho/data_spaces/phenotype_evaluator/phenotype_accumulator.hpp"

namespace clotho {
namespace genetics {

template < class BlockType, class WeightType, class TraitSpaceType >
class batch_phenotype_task< population_space< BlockType, WeightType >, TraitSpaceType > : public task {
public:
    typedef batch_phenotype_task< population_space< BlockType, WeightType >, TraitSpaceType >   self_type;

    typedef population_space< BlockType, WeightType >                                           space_type;
    typedef TraitSpaceType                                                                      trait_space_type;

    typedef typename space_type::genome_iterator                                                genome_iterator;

    typedef typename space_type::genome_type                                                haploid_genome_type;

    typedef typename space_type::base_genome_type::sequence_type::block_type                block_type;
    typedef typename space_type::base_genome_type::sequence_type::const_sequence_iterator   sequence_iterator;

    typedef typename trait_space_type::const_iterator                                       trait_iterator;

    typedef typename space_type::base_genome_type::sequence_type::bit_helper_type           bit_helper_type;

    typedef clotho::utility::debruijn_bit_walker< block_type >                              bit_walker_type;

    batch_phenotype_task( space_type * pop, trait_space_type * traits, unsigned int start, unsigned int end ) :
        m_pop( pop )
        , m_traits( traits )
        , m_start( start )
        , m_end( end )
    {}

    batch_phenotype_task( const self_type & other ) :
        m_pop( other.m_pop )
        , m_traits( other.m_traits )
        , m_start( other.m_start )
        , m_end( other.m_end )
    {}

    void operator()() {
        genome_iterator first = m_pop->begin_genomes(), last = m_pop->end_genomes();

        unsigned int i = 0;
        while( i < m_start ) { ++first; ++i; }

        while( i < m_end ) {
            if( first->first->isModifiable() ) {
                accumulate_allele_weights( first->first );
            }
            ++first;
            ++i;
        }
    }

    virtual ~batch_phenotype_task() {}

protected:

    void accumulate_allele_weights( const haploid_genome_type & genome ) {
        sequence_iterator first = genome->begin_sequence(), last = genome->end_sequence();

        typename space_type::base_genome_type::weight_type::weight_vector w( m_traits->trait_count(), 0 );

        size_t j = 0;
        while( first != last ) {
            block_type b = *first++;

            while( b ) {
                size_t bidx = j + bit_walker_type::unset_next_index( b );
                trait_iterator tfirst = m_traits->begin( bidx ), tlast = m_traits->end( bidx );

                unsigned int i = 0;
                while( tfirst != tlast ) {
                    w[ i++ ] += *tfirst++; 
                }
            }

            j += bit_helper_type::BITS_PER_BLOCK;
        }

        genome->update( w.begin(), w.end() );
    }

    space_type          * m_pop;
    trait_space_type    * m_traits;

    unsigned int        m_start, m_end;
};

}   // namespace genetics
}   // namespace clotho
#endif  // CLOTHO_BATCH_PHENOTYPE_TASK_POPULATION_SPACE_HPP_
