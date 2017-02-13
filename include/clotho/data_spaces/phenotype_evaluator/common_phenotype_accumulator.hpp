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
#ifndef CLOTHO_COMMON_PHENOTYPE_ACCUMULATOR_HPP_
#define CLOTHO_COMMON_PHENOTYPE_ACCUMULATOR_HPP_

#include "clotho/data_spaces/population_space/population_space_row.hpp"
#include "clotho/utility/debruijn_bit_walker.hpp"

namespace clotho {
namespace genetics {

template < class PopulationType, class TraitSpaceType >
class common_phenotype_accumulator;

template < class BlockType, class WeightType, class TraitSpaceType >
class common_phenotype_accumulator< population_space_row< BlockType, WeightType >, TraitSpaceType > {
public:
    typedef population_space_row< BlockType, WeightType >   space_type;

    typedef typename space_type::block_type                 block_type;

    typedef typename space_type::genome_pointer             genome_pointer;

    typedef TraitSpaceType                                  trait_space_type;

    typedef typename trait_space_type::weight_type          weight_type;

    typedef typename space_type::bit_helper_type                                            bit_helper_type;

    typedef clotho::utility::debruijn_bit_walker< block_type >                              bit_walker_type;

    common_phenotype_accumulator( trait_space_type * traits ) :
        m_traits( traits )
        , m_trait_count( traits->trait_count() )
        , m_weight_buffer( NULL )
    {
        m_weight_buffer = new weight_type[ m_trait_count ];
    }

    void operator()( genome_pointer first, genome_pointer last, block_type * neutrals ) {
        resetBuffer();

        block_type * f = neutrals;

        unsigned int j = 0;
        while( first != last ) {
            block_type n = *f++;
            block_type b = *first++;

            b &= (~n); 
            unsigned int bidx = j;
            while( b ) {
                bidx += bit_walker_type::next_and_shift(b);

                unsigned int k = 0, offset = bidx * m_trait_count;
                while( k < m_trait_count ) {
                    m_weight_buffer[ k++ ] += (*m_traits)[ offset++ ];
                }
            }

            j += bit_helper_type::BITS_PER_BLOCK;
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

#endif  // CLOTHO_COMMON_PHENOTYPE_ACCUMULATOR_HPP_
