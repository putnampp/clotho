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
#ifndef CLOTHO_POPULATION_SPACE_ROW_VECTOR_ARRAY_HPP_
#define CLOTHO_POPULATION_SPACE_ROW_VECTOR_ARRAY_HPP_

#include "clotho/data_spaces/phenotype_evaluator/trait_space_vector.hpp"

#include "clotho/utility/bit_helper.hpp"
#include <vector>

namespace clotho {
namespace genetics {

template < class BlockType, class WeightType >
class population_space_row {
public:
    typedef BlockType       block_type;

    typedef trait_space_vector< WeightType >           trait_space_type;

    typedef typename trait_space_type::weight_type      weight_type;
    typedef typename trait_space_type::vector_type      weight_vector;

    typedef typename weight_vector::iterator            weight_iterator;
    typedef typename weight_vector::const_iterator      const_weight_iterator;


    typedef clotho::utility::BitHelper< block_type >    bit_helper_type;

    typedef block_type *    individual_pointer;
    typedef block_type *    genome_pointer;

    typedef block_type *    row_pointer;

    population_space_row() :
        m_genetic_space( NULL )
        , m_genome_rows(0)
        , m_allele_block_columns(0)
        , m_allele_count(0)
        , m_trait_count(0)
        , m_max_rows(0)
        , m_max_cols(0)
    {}

    genome_pointer getHaploidGenome( unsigned int offset ) {
        return begin_genome( offset );
    }

//    individual_pointer getIndividual( unsigned int offset ) {
//        return getHaploidGenome( 2 * offset );
//    }

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

        m_allele_block_columns = (m_allele_count / bit_helper_type::BITS_PER_BLOCK) + ((m_allele_count % bit_helper_type::BITS_PER_BLOCK > 0) ? 1 : 0);
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

        assert( (m_genetic_space[ seq_idx ][ block_offset ] & mask) == (block_type)0);
#endif  // DEBUGGING

        m_genetic_space[ seq_idx ][ block_offset ] ^= mask;
    }

    void remove_fixed_allele( unsigned int all_idx ) {
#ifdef DEBUGGING
        assert( all_idx < m_allele_count );
#endif  // DEBUGGING

        unsigned int block_offset = (all_idx / bit_helper_type::BITS_PER_BLOCK );
        block_type mask = ~((block_type) 1 << (all_idx % bit_helper_type::BITS_PER_BLOCK ) );

        

        for( unsigned int i = 0; i < m_genome_rows; ++i ) {
#ifdef DEBUGGING
            assert( (m_genetic_space[ i ][ block_offset ] & ~mask) );
#endif  // DEBUGGING
            m_genetic_space[i][ block_offset ] &= mask;

        }
    }

    bool freeColumn(unsigned int idx ) {
        unsigned int block_offset = ( idx / bit_helper_type::BITS_PER_BLOCK );
        block_type mask = ((block_type) 1 << (idx % bit_helper_type::BITS_PER_BLOCK ));

        bool is_free = true;
        for( unsigned int i = 0; i < m_genome_rows; ++i ) {
            is_free = is_free && ((m_genetic_space[ i ][ block_offset ] & mask) == (block_type)0);
            if( !is_free ) {
                std::cerr << "Non-free space detected at index: " << idx << " in sequence " << i << std::endl;
                assert(false);
            }
        }

        return is_free;
    }

    row_pointer begin_block_row( unsigned int idx ) {
        return begin_genome(idx);
    }

    row_pointer end_block_row( unsigned int idx ) {
        return end_genome(idx);
    }

    genome_pointer begin_genome( unsigned int idx ) {
        assert( idx < m_genome_rows );
        return m_genetic_space[ idx ];
    }

    genome_pointer end_genome( unsigned int idx ) {
        assert( idx < m_genome_rows );

        return m_genetic_space[ idx ] + m_allele_block_columns;
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

    const_weight_iterator begin_genome_traits( unsigned int idx ) const {
        return m_genome_traits.begin() + idx * m_trait_count;
    }

    const_weight_iterator end_genome_traits( unsigned int idx ) const {
        return m_genome_traits.begin() + (idx + 1) * m_trait_count;
    }

    virtual ~population_space_row() {
        if( m_genetic_space != NULL ) {
            for( unsigned int i = 0; i < m_max_rows; ++i )
                delete [] m_genetic_space[ i ];

            delete [] m_genetic_space;
        }
    }

    void clear() {
        for( unsigned int i = 0; i < m_max_rows; ++i )
            memset( m_genetic_space[i], 0, sizeof(block_type) * m_allele_block_columns );
    }

protected:

    void resize( ) {

        if( m_allele_block_columns > m_max_cols ) {
            // if genome_rows == max_rows && block_columns > max_cols then re-allocate all rows
            // if genome_rows < max_rows && block_columns > max_cols then re-allocate all rows
            // if genome_rows > max_rows && block_columns > max_cols then re-allocate all rows
            if( m_genetic_space != NULL ) {
                for( unsigned int i = 0; i < m_max_rows; ++i )
                    delete [] m_genetic_space[ i ];

                delete [] m_genetic_space;
            }

            m_max_cols = m_allele_block_columns + 100;  // add 100 blocks/row in an attempt to reduce malloc calls (6400 bits)

            m_max_rows = m_genome_rows;

            m_genetic_space = new block_type*[ m_max_rows ];
            for( unsigned int i = 0; i < m_max_rows; ++i )
                m_genetic_space[ i ] = new block_type[ m_max_cols ];

        } else if( m_genome_rows > m_max_rows ) {
            // if genome_rows > max_rows && block_columns == max_cols then append new rows ( copy and extend )
            // if genome_rows > max_rows && block_columns < max_cols then append new rows ( copy and extend using current max_cols )
            //
            block_type ** tmp = new block_type*[ m_genome_rows ];

            for( unsigned int i = 0; i < m_max_rows; ++i )
                tmp[ i ] = m_genetic_space[ i ];

            for( unsigned int i = m_max_rows; i < m_genome_rows; ++i )
                tmp[ i ] = new block_type[ m_max_cols ];

            if( m_genetic_space != NULL ) {
                delete [] m_genetic_space;
            }

            m_genetic_space = tmp;
            m_max_rows = m_genome_rows;
        }
        // else {
        // if genome_rows == max_rows && block_columns == max_cols then do nothing
        // if genome_rows < max_rows && block_columns == max_cols then do nothing ( leave existing max rows allocated )
        // if genome_rows == max_rows && block_columns < max_cols then do nothing ( leave existing max columns allocated )
        // if genome_rows < max_rows && block_columns < max_cols then do nothing ( leave existing maxs allocated )
        // }

        unsigned int new_total = m_trait_count * m_max_rows;

        while( m_genome_traits.size() < new_total ) {
            m_genome_traits.push_back(0);
        }
    }

    block_type      ** m_genetic_space;

    unsigned int    m_genome_rows, m_allele_block_columns, m_allele_count, m_trait_count;

    unsigned int    m_max_rows, m_max_cols;

    weight_vector   m_genome_traits;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_POPULATION_SPACE_ROW_VECTOR_ARRAY_HPP_
