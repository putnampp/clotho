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
#ifndef CLOTHO_COMMON_PHENOTYPE_ACCUMULATOR_INDEX_VECTOR_ALIGNMENT_HPP_
#define CLOTHO_COMMON_PHENOTYPE_ACCUMULATOR_INDEX_VECTOR_ALIGNMENT_HPP_

#include "clotho/data_spaces/population_space/population_spaces.hpp"
#include "clotho/data_spaces/phenotype_evaluator/common_phenotype_accumulator_def.hpp"
#include "clotho/utility/debruijn_bit_walker.hpp"

namespace clotho {
namespace genetics {

template < class IndexType, class BlockType, class WeightType >
class common_phenotype_accumulator< population_space< index_vector_alignment< IndexType, BlockType >, trait_space_vector< WeightType > >, trait_space_vector< WeightType > > {
public:
    typedef population_space< index_vector_alignment< IndexType, BlockType >, trait_space_vector< WeightType > >   space_type;
    typedef trait_space_vector< WeightType >                trait_space_type;

    typedef typename space_type::index_type                 index_type;
    typedef typename space_type::genome_pointer             genome_pointer;

    typedef typename trait_space_type::weight_type          weight_type;

    common_phenotype_accumulator( trait_space_type * traits ) :
        m_traits( traits )
        , m_trait_count( traits->trait_count() )
        , m_weight_buffer( NULL )
    {
        m_weight_buffer = new weight_type[ m_trait_count ];
    }

    template< class AlleleType >
    void operator()( genome_pointer first, genome_pointer last, AlleleType * all ) {
        resetBuffer();

        while( first != last ) {
            index_type idx = *first++;

            if( !all->isNeutral( idx ) ) {

                unsigned int k = 0, offset = idx * m_trait_count;
                while( k < m_trait_count ) {
                    m_weight_buffer[ k++ ] = (*m_traits)[ offset++ ];
                }
            }
        }
    }

    void resetBuffer( ) {
        for( unsigned int i = 0; i < m_trait_count; ++i ) {
            m_weight_buffer[ i ] = 0;
        }
    }

    weight_type * getResults() {
        return m_weight_buffer;
    }

    unsigned int trait_count() const {
        return m_trait_count;
    }

    virtual ~common_phenotype_accumulator() {
        if( m_weight_buffer != NULL ) {
            delete [] m_weight_buffer;
        }
    }

protected:
    trait_space_type * m_traits;
    unsigned int m_trait_count;
    weight_type * m_weight_buffer;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_COMMON_PHENOTYPE_ACCUMULATOR_INDEX_VECTOR_ALIGNMENT_HPP_
