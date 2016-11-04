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
#ifndef CLOTHO_BATCH_PHENOTYPE_TASK_POPULATION_SPACE_ROW_HPP_
#define CLOTHO_BATCH_PHENOTYPE_TASK_POPULATION_SPACE_ROW_HPP_

#include "clotho/data_spaces/task/task.hpp"

#include "clotho/utility/debruijn_bit_walker.hpp"
#include "clotho/data_spaces/phenotype_evaluator/phenotype_accumulator.hpp"

#include "clotho/data_spaces/population_space/population_space_row.hpp"

namespace clotho {
namespace genetics {

template < class BlockType, class WeightType, class TraitSpaceType >
class batch_phenotype_task< population_space_row< BlockType, WeightType >, TraitSpaceType > : public task {
public:
    typedef batch_phenotype_task< population_space_row< BlockType, WeightType >, TraitSpaceType >   self_type;

    typedef population_space_row< BlockType, WeightType >                                           space_type;
    typedef TraitSpaceType                                                                      trait_space_type;

    typedef typename space_type::genome_pointer                                             genome_pointer;

    typedef typename space_type::block_type                                                 block_type;

    typedef typename space_type::weight_vector                                              weight_vector;

    typedef typename trait_space_type::const_iterator                                       trait_iterator;

    typedef typename space_type::bit_helper_type                                            bit_helper_type;

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
        const unsigned int TRAIT_COUNT = m_traits->trait_count();

        weight_vector local_weights( TRAIT_COUNT, 0 );
        for( unsigned int i = m_start; i < m_end; ++i ) {
            genome_pointer first = m_pop->begin_genome(i), last = m_pop->end_genome(i);

            accumulate_allele_weights_alt2( first, last, local_weights );

            m_pop->updateGenomeWeights( i, local_weights );

            for( unsigned int j = 0; j < TRAIT_COUNT; ++j ) {
                local_weights[ j ] = 0.0;
            }
        }
    }

    virtual ~batch_phenotype_task() {}

protected:
    
    void accumulate_allele_weights( genome_pointer first, genome_pointer last, weight_vector & weights ) {
        unsigned int j = 0;
        while( first != last ) {
            block_type b = *first++;

            unsigned int bidx = j;
            while( b ) {
                bidx += bit_walker_type::next_and_shift(b);

                trait_iterator tfirst = m_traits->begin(bidx), tlast = m_traits->end(bidx);

                unsigned int k = 0;
                while( tfirst != tlast ) {
                    weights[ k++ ] += *tfirst++;
                }
            }

            j += bit_helper_type::BITS_PER_BLOCK;
        }
    }

    void accumulate_allele_weights_alt( genome_pointer first, genome_pointer last, weight_vector & weights ) {
        const unsigned int TRAIT_COUNT = m_traits->trait_count();
        unsigned int j = 0;

        trait_iterator tfirst = m_traits->begin(0);

        while( first != last ) {
            block_type b = *first++;

            unsigned int bidx = j;
            while( b ) {
                bidx += bit_walker_type::next_and_shift(b);

                trait_iterator tmp = (tfirst + bidx * TRAIT_COUNT);
                unsigned int k = 0;
                while( k < TRAIT_COUNT ) {
                    weights[ k++ ] += *tmp++;
                }
            }

            j += bit_helper_type::BITS_PER_BLOCK;
        }
    }

    void accumulate_allele_weights_alt2( genome_pointer first, genome_pointer last, weight_vector & weights ) {
        const unsigned int TRAIT_COUNT = m_traits->trait_count();
        unsigned int j = 0;

        unsigned int index_buffer[ bit_helper_type::BITS_PER_BLOCK ];

        while( first != last ) {
            block_type b = *first++;

            unsigned int offset = 0;

            unsigned int bidx = j;
            while( b ) {
                bidx += bit_walker_type::next_and_shift(b);
                index_buffer[ offset++ ] = bidx;
            }

            while( offset-- ) {
                unsigned int k = 0, l = index_buffer[ offset ] * TRAIT_COUNT;
                while( k < TRAIT_COUNT ) {
                    weights[ k++ ] += (*m_traits)[ l++ ];
                }
            }

            j += bit_helper_type::BITS_PER_BLOCK;
        }
    }

    space_type          * m_pop;
    trait_space_type    * m_traits;

    unsigned int        m_start, m_end;
};

}   // namespace genetics
}   // namespace clotho
#endif  // CLOTHO_BATCH_PHENOTYPE_TASK_POPULATION_SPACE_ROW_HPP_
