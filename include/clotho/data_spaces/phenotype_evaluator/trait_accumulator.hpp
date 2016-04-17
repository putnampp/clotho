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
#ifndef CLOTHO_TRAIT_WEIGHT_ACCUMULATOR_HPP_
#define CLOTHO_TRAIT_WEIGHT_ACCUMULATOR_HPP_

#include <vector>

#include "clotho/data_spaces/trait_type_of.hpp"
#include "clotho/utility/debruijn_bit_walker.hpp"

template < class GeneticSpaceType >
class TraitWeightAccumulator {
public:
    
    typedef GeneticSpaceType    genetic_space_type;

    typedef typename genetic_space_type::allele_type            allele_type;
    
    typedef clotho::genetics::trait_type_of< allele_type >      trait_helper_type;
    typedef typename trait_helper_type::type                    trait_vector_type;

    typedef std::vector< trait_vector_type >    accumulator_type;

    typedef clotho::utility::debruijn_bit_walker< block_type >   bit_walker_type;

    TraitWeightAccumulator( const genetic_space_type & genomes ) {
        update( genomes );
    }

    void update( const genetic_space_type & genomes ) {

        m_trait_weights.reserve( genomes.sequence_count() );

        typedef typename genetic_space_type::const_block_iterator   const_block_iterator;
        typedef typename genetic_space_type::block_type             block_type;

        typedef typename allele_type::const_trait_iterator   const_trait_iterator;
        typedef typename allele_type::weight_type            weight_type;

        const_block_iterator block_iter = genomes.block_iterator();
        size_t i = 0, j = 0;
        while( block_iter.hasNext() ) {
            block_type b = block_iter.next();

            if( j == 0 ) {
                if( j >= m_trait_weights.size() ) {
                    m_trait_weights.push_back( trait_helper_type::makeEmptyTrait( genomes ) );
                } else {
                    trait_helper_type::resetTrait(m_trait_weights[j])
                }
            }

            trait_type & t = m_trait_weights[ j ];

            while( b ) {
                unsigned int b_idx = bit_walker_type::unset_next_index( b );
                const_trait_iterator trait_it = genomes.getTraitIterator( j + b_idx );
                size_t k = 0;
                while( trait_it.hasNext() ) {
                    weight_type w = trait_it.next();
                    t[ k++ ] += w;
                }
            }

            if( ++i >= N ) {
                i = 0;
                j += genetic_space_type::bit_helper_type::BITS_PER_BLOCK;
            }
        }
    }

    trait_type &  getTraitAt( size_t idx ) {
        assert( 0 <= idx && idx < m_trait_weights.size() );
        return m_trait_weights[ idx ];
    }

    virtual ~TraitWeightAccumulator() {}
protected:

    accumulator_type   m_trait_weights;    
};

#endif  // CLOTHO_TRAIT_WEIGHT_ACCUMULATOR_HPP_
