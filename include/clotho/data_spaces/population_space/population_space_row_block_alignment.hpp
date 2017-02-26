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
#ifndef CLOTHO_POPULATION_SPACE_ROW_BLOCK_ALIGNMENT_HPP_
#define CLOTHO_POPULATION_SPACE_ROW_BLOCK_ALIGNMENT_HPP_

#include "clotho/data_spaces/population_space/population_space_def.hpp"
#include "clotho/data_spaces/phenotype_evaluator/trait_space_vector.hpp"

#include "clotho/utility/bit_helper.hpp"
#include <vector>

namespace clotho {
namespace genetics {

template < class BlockType >
struct row_block_alignment {
    typedef BlockType   block_type;
};

template < class BlockType, class WeightType >
class population_space< row_block_alignment< BlockType >, trait_space_vector< WeightType > > {
public:
    typedef row_block_alignment< BlockType >    alignment_type;
    typedef trait_space_vector< WeightType >    trait_space_type;

    typedef typename alignment_type::block_type block_type;

    typedef typename trait_space_type::weight_type      weight_type;
    typedef typename trait_space_type::vector_type      weight_vector;

    typedef typename weight_vector::iterator            weight_iterator;
    typedef typename weight_vector::const_iterator      const_weight_iterator;

    typedef clotho::utility::BitHelper< block_type >    bit_helper_type;

    typedef block_type *    individual_pointer;
    typedef block_type *    genome_pointer;

    typedef block_type *    row_pointer;

    population_space() :
        m_genetic_space( NULL )
        , m_genetic_space_allocated(0)
        , m_genome_rows(0)
        , m_allele_block_columns(0)
        , m_allele_count(0)
        , m_trait_count(0)
    {}

    genome_pointer getHaploidGenome( unsigned int offset ) {
        assert( offset < m_genome_rows );
        return m_genetic_space + offset * m_allele_block_columns;;
    }

    individual_pointer getIndividual( unsigned int offset ) {
        return getHaploidGenome( 2 * offset );
    }

    void grow( unsigned int I, unsigned int A, unsigned int T ) {
        setIndividualCount( I );
        setMaxAlleles( A );
        setTraitCount( T );

        resize( );

        //clear();
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

//        m_allele_block_columns = (m_allele_count / bit_helper_type::BITS_PER_BLOCK) + ((m_allele_count % bit_helper_type::BITS_PER_BLOCK > 0) ? 1 : 0);
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
        unsigned int block_offset = (all_idx / bit_helper_type::BITS_PER_BLOCK);
        unsigned int bit_offset= (all_idx % bit_helper_type::BITS_PER_BLOCK);

        block_type mask = ((block_type)1 << bit_offset);

#ifdef DEBUGGING
        assert( seq_idx < m_genome_rows );
        assert( all_idx < m_allele_count);

        assert( (m_genetic_space[ block_offset + seq_idx * m_allele_block_columns ] & mask) == (block_type)0);
#endif  // DEBUGGING

        m_genetic_space[ block_offset + seq_idx * m_allele_block_columns ] ^= mask;
    }

    void remove_fixed_allele( unsigned int all_idx ) {
#ifdef DEBUGGING
        assert( all_idx < m_allele_count );
#endif  // DEBUGGING

        unsigned int block_offset = (all_idx / bit_helper_type::BITS_PER_BLOCK );
        block_type mask = ~((block_type) 1 << (all_idx % bit_helper_type::BITS_PER_BLOCK ) );

        const unsigned int STEP = getMaxBlocks();
        
        block_type * first = m_genetic_space + block_offset;

        for( unsigned int i = 0; i < m_genome_rows; ++i ) {
#ifdef DEBUGGING
            assert( (*first & ~mask) );
#endif  // DEBUGGING
            *first &= mask;

            first += STEP;
        }
    }

    bool freeColumn(unsigned int idx ) {
        unsigned int block_offset = ( idx / bit_helper_type::BITS_PER_BLOCK );
        block_type mask = ((block_type) 1 << (idx % bit_helper_type::BITS_PER_BLOCK ));

        block_type * first = m_genetic_space + block_offset;

        bool is_free = true;
        for( unsigned int i = 0; i < m_genome_rows; ++i ) {
            is_free = is_free && ((*first & mask) == (block_type)0);
            if( !is_free ) {
                std::cerr << "Non-free space detected at index: " << idx << " in sequence " << i << std::endl;
                assert(false);
            }
        }

        return is_free;
    }

    row_pointer begin_block_row( unsigned int idx ) {
        return m_genetic_space + idx * m_allele_block_columns;
    }

    row_pointer end_block_row( unsigned int idx ) {

        return m_genetic_space + (idx + 1) * m_allele_block_columns;
    }

    genome_pointer begin_genome( unsigned int idx ) {
        if( idx >= m_genome_rows )
            return m_genetic_space;
    
        return m_genetic_space + idx * m_allele_block_columns;
    }

    genome_pointer end_genome( unsigned int idx ) {
        if( idx >= m_genome_rows )
            return m_genetic_space;

        return m_genetic_space + (idx + 1) * m_allele_block_columns;
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
        memset( m_genetic_space, 0, sizeof(block_type) * m_genetic_space_allocated );
    }

    virtual ~population_space() {
        if( m_genetic_space != NULL ) {
            delete [] m_genetic_space;
        }
    }

protected:

    void resize( ) {

        size_t  new_total = m_genome_rows * m_allele_block_columns;

        if( new_total > m_genetic_space_allocated ) {
            if( m_genetic_space != NULL ) {
                delete [] m_genetic_space;
            }

            //new_total += 10 * m_genome_rows;  // pad each allocation by an additional 10 rows
            m_genetic_space = new block_type[ new_total ];

            m_genetic_space_allocated = new_total;
        }

        new_total = m_trait_count * m_genome_rows;

        while( m_genome_traits.size() < new_total ) {
            m_genome_traits.push_back(0);
        }
    }

    block_type      * m_genetic_space;

    size_t          m_genetic_space_allocated;

    unsigned int    m_genome_rows, m_allele_block_columns, m_allele_count, m_trait_count;

    weight_vector   m_genome_traits;
};

}   // namespace genetics
}   // namespace clotho

namespace clotho {
namespace utility {

template < class BlockType, class WeightType >
struct state_getter< clotho::genetics::population_space< clotho::genetics::row_block_alignment< BlockType >, clotho::genetics::trait_space_vector< WeightType > > > {
    typedef clotho::genetics::population_space< clotho::genetics::row_block_alignment< BlockType >, clotho::genetics::trait_space_vector< WeightType > > object_type;

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

#endif  // CLOTHO_POPULATION_SPACE_ROW_BLOCK_ALIGNMENT_HPP_

