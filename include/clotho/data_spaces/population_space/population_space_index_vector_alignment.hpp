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
#ifndef CLOTHO_POPULATION_SPACE_VECTOR_ALIGNMENT_HPP_
#define CLOTHO_POPULATION_SPACE_VECTOR_ALIGNMENT_HPP_

#include "clotho/data_spaces/population_space/population_space_def.hpp"
#include "clotho/data_spaces/phenotype_evaluator/trait_space_vector.hpp"

#include "clotho/utility/bit_helper.hpp"

#include <vector>
#include <algorithm>

namespace clotho {
namespace genetics {

template < class IndexType, class BlockType >
struct index_vector_alignment {
    typedef IndexType   index_type;
    typedef BlockType   block_type;
};

template < class IndexType, class BlockType, class WeightType >
class population_space< index_vector_alignment< IndexType, BlockType >, trait_space_vector< WeightType > > {
public:
    typedef index_vector_alignment< IndexType, BlockType >         alignment_type;
    typedef trait_space_vector< WeightType >            trait_space_type;

    typedef typename alignment_type::index_type         index_type;
    typedef typename alignment_type::block_type         block_type;

    typedef typename trait_space_type::weight_type      weight_type;
    typedef typename trait_space_type::vector_type      weight_vector;

    typedef typename weight_vector::iterator            weight_iterator;
    typedef typename weight_vector::const_iterator      const_weight_iterator;

    typedef clotho::utility::BitHelper< block_type >    bit_helper_type;

    typedef std::vector< index_type >                   index_vector_type;
    typedef std::vector< index_vector_type >            index_space_type;

    typedef typename index_vector_type::iterator        individual_pointer;
    typedef typename index_vector_type::iterator        genome_pointer;

    typedef typename index_space_type::iterator         row_pointer;

    population_space() :
        m_genetic_space_allocated(0)
        , m_genome_rows(0)
        , m_allele_block_columns(0)
        , m_allele_count(0)
        , m_trait_count(0)
    {}

    genome_pointer getHaploidGenome( unsigned int offset ) {
        return begin_genome( offset );
    }

    individual_pointer getIndividual( unsigned int offset ) {
        return getHaploidGenome( 2 * offset );
    }

    void grow( unsigned int I, unsigned int A, unsigned int T ) {
        setIndividualCount( I );
        setMaxAlleles( A );
        setTraitCount( T );

        resize( );
    }

    unsigned int getIndividualCount( ) const {
        return (m_genome_rows / 2);
    }

    void setIndividualCount( unsigned int I ) {
        m_genome_rows = 2 * I;
    }

    unsigned int haploid_genome_count() const {
        return m_genome_rows;
    }

    unsigned int getMaxBlocks() const {
        return m_allele_block_columns;
    }

    unsigned int getMaxAlleles( ) const {
        return m_allele_count;
    }

    void setMaxAlleles( unsigned int A ) {
        m_allele_count = A;

        unsigned int allele_cache_lines = (m_allele_count / bit_helper_type::BITS_PER_CACHE_LINE) + ((m_allele_count % bit_helper_type::BITS_PER_CACHE_LINE > 0) ? 1 : 0);

        m_allele_block_columns = allele_cache_lines * bit_helper_type::BLOCKS_PER_CACHE_LINE;
    }

    unsigned int getTraitCount() const {
        return m_trait_count;
    }

    void setTraitCount( unsigned int T ) {
        m_trait_count = T;
    }

    void mutate( unsigned int seq_idx, unsigned int all_idx ) {
        assert( seq_idx < m_genome_rows );
        assert( all_idx < m_allele_count);
        
        m_genetic_space[ seq_idx ].push_back( all_idx );

        std::sort( m_genetic_space[ seq_idx ].begin(), m_genetic_space[ seq_idx ].end() );
    }

    void remove_fixed_allele( unsigned int all_idx ) {

        unsigned int N = 0;
        for( row_pointer it = m_genetic_space.begin(); it != m_genetic_space.end(); it++ ) {
            genome_pointer git = std::lower_bound( it->begin(), it->end(), all_idx );

            if( git != it->end() ) {
                it->erase( git );
                ++N;
            }
        }

        assert( N == m_genetic_space.size() );
    }

    bool freeColumn(unsigned int idx ) {
        bool is_free = true;
        for( row_pointer it = m_genetic_space.begin(); it != m_genetic_space.end() && is_free; it++ ) {
            genome_pointer git = std::lower_bound( it->begin(), it->end(), idx );

            is_free = (git == it->end());
        }

        return is_free;
    }

    row_pointer begin_row( unsigned int idx ) {
        assert( idx < m_genetic_space.size() );

        return m_genetic_space.begin() + idx;
    }

    row_pointer end_row( unsigned int idx ) {
        if( idx + 1 >= m_genetic_space.size() ) {
            return m_genetic_space.end();
        }
        return m_genetic_space.begin() + (idx + 1);
    }

    genome_pointer begin_genome( unsigned int idx ) {
        assert( idx < m_genetic_space.size() );
    
        return m_genetic_space[ idx ].begin();
    }

    genome_pointer end_genome( unsigned int idx ) {
        assert( idx < m_genetic_space.size() );

        return m_genetic_space[ idx ].end();
    }

    void updateGenomeWeights( unsigned int genome_index, weight_vector & weights ) {
#ifdef DEBUGGING
        assert( genome_index < m_genome_rows );
        assert( weights.size() <= m_trait_count );
#endif  // DEBUGGING

        const_weight_iterator w = weights.begin();

        weight_iterator first = m_genome_traits.begin() + genome_index * m_trait_count;

        while( w != weights.end() ) {
            *first++ = *w++;
        }
    }

    void updateGenomeWeights( unsigned int genome_index, weight_type * weights ) {
        for(unsigned int i = 0, j = genome_index * m_trait_count; i < m_trait_count; ++i, ++j ) {
            m_genome_traits[ j ] = weights[ i ];
        }
    }

    const_weight_iterator begin_genome_traits( unsigned int idx ) const {
        return m_genome_traits.begin() + idx * m_trait_count;
    }

    const_weight_iterator end_genome_traits( unsigned int idx ) const {
        return m_genome_traits.begin() + (idx + 1) * m_trait_count;
    }

    void clear() {
        for( row_pointer it = m_genetic_space.begin(); it != m_genetic_space.end(); it++ ) {
            it->clear();
        }
    }

    virtual ~population_space() {
//        if( m_genetic_space != NULL ) {
//            delete [] m_genetic_space;
//        }
    }

protected:

    void resize( ) {

        while( m_genetic_space.size() < m_genome_rows ) {
            m_genetic_space.push_back( index_vector_type() );
        }

        size_t new_total = m_trait_count * m_genome_rows;

        while( m_genome_traits.size() < new_total ) {
            m_genome_traits.push_back(0);
        }
    }

    index_space_type m_genetic_space;

    size_t          m_genetic_space_allocated;

    unsigned int    m_genome_rows, m_allele_block_columns, m_allele_count, m_trait_count;

    weight_vector   m_genome_traits;
};

}   // namespace genetics
}   // namespace clotho

namespace clotho {
namespace utility {

template < class IndexType, class BlockType, class WeightType >
struct state_getter< clotho::genetics::population_space< clotho::genetics::index_vector_alignment< IndexType, BlockType >, clotho::genetics::trait_space_vector< WeightType > > > {
    typedef clotho::genetics::population_space< clotho::genetics::index_vector_alignment< IndexType, BlockType >, clotho::genetics::trait_space_vector< WeightType > > object_type;

    void operator()( boost::property_tree::ptree & s, object_type & obj ) {
        boost::property_tree::ptree ph;

        for( unsigned int i = 0; i < obj.haploid_genome_count(); ++i ) {
            boost::property_tree::ptree p;

            clotho::utility::add_value_array( p, obj.begin_genome_traits( i ), obj.end_genome_traits( i ) );

            ph.push_back( std::make_pair( "", p ) );
        }

        s.put_child( "phenotype", ph );
    }
};

}   // namespace utility
}   // namespace clotho

#endif  // CLOTHO_POPULATION_SPACE_VECTOR_ALIGNMENT_HPP_

