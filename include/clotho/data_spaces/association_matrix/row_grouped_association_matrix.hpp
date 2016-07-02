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
#ifndef CLOTHO_ROW_GROUPED_ASSOCIATION_MATRIX_HPP_
#define CLOTHO_ROW_GROUPED_ASSOCIATION_MATRIX_HPP_

#include "clotho/data_spaces/association_matrix/association_matrix_def.hpp"
#include "clotho/data_spaces/association_matrix/types/row_grouped.hpp"

#include "clotho/data_spaces/growable2D.hpp"
#include "clotho/utility/bit_helper.hpp"
#include "clotho/data_spaces/association_matrix/association_vector.hpp"

#ifdef DEBUGGING
#define ASSERT_VALID_RANGE( x, min, max) assert( min <= x && x < max);
#else   // DEBUGGING
#define ASSERT_VALID_RANGE( x, min, max)
#endif  // DEBUGGING

namespace clotho {
namespace genetics {

template < class BlockType, unsigned char P >
class association_matrix< BlockType, row_grouped< P > > : public growable2D {
public:
    typedef BlockType                                   block_type;
    typedef row_grouped< P >                            alignment_type;

    typedef clotho::utility::BitHelper< block_type >    bit_helper_type;

    typedef block_type *                                raw_block_pointer;

    association_matrix( size_t rows = 1, size_t columns = 1 ) : 
        m_data( NULL )
        , m_rows( 0 )
        , m_columns( 0 )
        , m_allocated_size( 0 )
        , m_bpr(0)
    {
        this->resize( rows, columns );
    }

    size_t row_count( ) const {
        return m_rows;
    }

    size_t column_count() const {
        return m_columns;
    }

    bool operator()( size_t row, size_t col ) {
        ASSERT_VALID_RANGE( row, 0, m_rows );
        ASSERT_VALID_RANGE( col, 0, m_columns );

        size_t idx = block_index( row, col, m_bpr );
        return ((m_data[ idx ] & bit_helper_type::bit_offset( col )) != 0);
    }

    void flip( size_t row, size_t col ) {
        ASSERT_VALID_RANGE( row, 0, m_rows );
        ASSERT_VALID_RANGE( col, 0, m_columns );

        size_t idx = block_index( row, col, m_bpr );
#ifdef DEBUGGING
        assert( (m_data[ idx ] & bit_helper_type::bit_offset( col )) == 0 );
#endif  // DEBUGGING

        m_data[ idx ] ^= bit_helper_type::bit_offset( col );
    }

    bool freeColumn( size_t fr_idx ) {
        block_type mask = bit_helper_type::bit_offset( fr_idx );

        size_t i = 0;
        bool is_free = true;
        while( is_free && i < m_rows ) {
            size_t idx = block_index( i++, fr_idx, m_bpr );

            is_free = ((m_data[idx] & mask) == 0);
            if( !is_free) std::cerr << "Row " << (i - 1) << " at " << idx << " [" << (m_data + idx) << "] -> " << std::hex << (m_data[idx] & mask) << std::dec << std::endl;
        }
        return is_free;
    }

    void flipColumn( size_t fixed_index ) {
        block_type mask = bit_helper_type::bit_offset( fixed_index );

        size_t i = 0;
        while( i < m_rows ) {
            size_t idx = block_index( i++, fixed_index, m_bpr );

#ifdef DEBUGGING
            assert( (m_data[idx] & mask) != 0 );
#endif  // DEBUGGING
            m_data[ idx ] ^= mask;
        }
    }

    size_t block_row_count() const {
        return m_rows;
    }

    size_t block_column_count() const {
        return m_bpr;
    }

    size_t size() const {
        return m_columns * m_rows;
    }

    size_t allocated_size() const {
        return m_allocated_size;
    }

    raw_block_pointer begin_row( size_t idx ) const {
        ASSERT_VALID_RANGE( idx, 0, row_count() );
        size_t offset = block_row_offset( idx, m_bpr );
        return m_data + offset;
    }

    raw_block_pointer end_row( size_t idx ) const {
        ASSERT_VALID_RANGE( idx, 0, row_count() );
        return m_data + block_row_offset( idx + alignment_type::GROUP_SIZE, m_bpr );
    }

    raw_block_pointer begin_block_row( size_t idx ) const {
        ASSERT_VALID_RANGE( idx, 0, block_row_count() );

        return m_data + block_row_offset( idx, m_bpr);
    }

    raw_block_pointer end_block_row( size_t idx ) const {
        ASSERT_VALID_RANGE( idx, 0, block_row_count() );
        return m_data + block_row_offset( idx + alignment_type::GROUP_SIZE, m_bpr );
    }

    virtual size_t grow( size_t rows, size_t cols ) {
        this->resize( rows, cols );
        return this->size();
    }

    void clear( ) {
        memset( m_data, 0, sizeof( block_type ) * m_allocated_size );
    }

    virtual ~association_matrix() {
        if( m_data != NULL ) {
            delete [] m_data;
        }
    }

protected:

/**
 * bpr - blocks_per_row
 */
    inline size_t block_index( size_t row, size_t col, size_t bpr ) const {
        return block_row_offset( row, bpr ) + block_column_offset( col );
    }

    inline size_t block_row_offset( size_t row, size_t bpr ) const {
        return (row / row_grouped< P >::GROUP_SIZE) * bpr * row_grouped< P >::GROUP_SIZE + (row % row_grouped< P >::GROUP_SIZE);
    }

    inline size_t block_column_offset( size_t col ) const {
        return (col / clotho::utility::BitHelper< BlockType >::BITS_PER_BLOCK) * row_grouped< P >::GROUP_SIZE;
    }

    void resize( size_t r, size_t c ) {
        size_t blocks_per_row =  bit_helper_type::padded_block_count( c );

        size_t padded_rows = r + (row_grouped< P >::GROUP_SIZE - (r % row_grouped< P >::GROUP_SIZE) );
        size_t new_size = padded_rows * blocks_per_row;

        if( new_size > m_allocated_size ) {
            if( m_data != NULL ) {
                delete [] m_data;
            }

            m_data = new block_type[ new_size ];

            m_allocated_size = new_size;

            assert( m_data != NULL );

            clear();
        }

        m_rows = r;
        m_columns = c;
        m_bpr = blocks_per_row;
    }

    raw_block_pointer   m_data;

    size_t m_rows, m_columns, m_allocated_size;
    size_t m_bpr;
};

template < class BlockType >
class association_matrix< BlockType, row_grouped< 1 > > : public growable2D {
public:
    typedef BlockType                                   block_type;
    typedef row_grouped< 1 >                            alignment_type;

    typedef clotho::utility::BitHelper< block_type >    bit_helper_type;

    typedef association_vector< BlockType >             row_vector;
    typedef std::vector< row_vector >                   data_vector;

    typedef typename data_vector::iterator              row_iterator;

    typedef typename row_vector::const_raw_pointer      raw_block_pointer;

    association_matrix( size_t rows = 1, size_t columns = 1 ) : 
//        m_data( NULL )
         m_rows( 0 )
        , m_columns( 0 )
        , m_bpr(0)
    {
        this->resize( rows, columns );
    }

    size_t row_count( ) const {
        return m_rows;
    }

    size_t column_count() const {
        return m_columns;
    }

    bool operator()( size_t row, size_t col ) {
        ASSERT_VALID_RANGE( row, 0, m_rows )
        ASSERT_VALID_RANGE( col, 0, m_columns )

        size_t idx = block_index( row, col, m_bpr );
        return ((m_data[ idx ] & bit_helper_type::bit_offset( col )) != 0);
    }

    void flip( size_t row, size_t col ) {
        ASSERT_VALID_RANGE( row, 0, m_rows )
        ASSERT_VALID_RANGE( col, 0, m_columns )

        if( m_data[ row ].m_readonly ) {
            row_vector a( m_bpr );

            a.copy( m_data[row] );

            m_data[ row ] = a;
        }

        size_t offset = block_column_offset(col);
#ifdef DEBUGGING
        assert( !m_data[ row ].m_readonly );
        assert( offset < m_data[ row ].m_size );

        if( col == 5593 ) {
            std::cerr << "Setting column: " << col << "; row: " << row << std::endl;
        }
#endif  // DEBUGGING

        m_data[ row ].get()[ offset ]  ^= bit_helper_type::bit_offset( col );
    }

#ifdef DEBUGGING
    bool freeColumn( size_t fr_idx ) {
        block_type mask = bit_helper_type::bit_offset( fr_idx );
        size_t offset = block_column_offset( fr_idx );
        size_t i = 0;
        bool is_free = true;
        while( is_free && i < m_rows ) {
            if( offset < m_data[ i ].m_size ) {
                block_type b = m_data[ i ].get()[ offset ];
                is_free = ((b & mask) == 0);
                if( !is_free)
                    std::cerr << "Row " << i << " at " << offset << " [" << m_data[i].get() + offset << "] -> is set" << std::endl;
            }
            ++i;
        }
        return is_free;
    }
#endif  // DEBUGGING

    void flipColumn( size_t fixed_index ) {
        block_type mask = bit_helper_type::bit_offset( fixed_index );
        size_t offset = block_column_offset( fixed_index );

        for( row_iterator it = m_data.begin(); it != m_data.end(); ++it ) {
#ifdef DEBUGGING
            assert( offset < it->m_size );
#endif  // DEBUGGING
            block_type b = it->get()[ offset ];
            if( (b & mask) != 0 ) {
                b ^= mask;
                it->get()[ offset ] = b;
            }
        }
    }

    void finalize() {
#ifdef DEBUGGING
        assert( m_data.size() == m_rows );
#endif  // DEBUGGING
        for( size_t i = 0; i < m_data.size(); ++i )
            m_data[i].finalize();
    }

    void fill_empty() {
        row_vector a( m_bpr ); // create an empty row vector of desired size

        a.finalize();

        for( size_t i = m_data.size(); i < m_rows; ++i ) {
            m_data.push_back( row_vector( a ) ); // fill the space with copies of the empty vector

        }
    }

    void push_back( row_vector & r ) {
        m_data.push_back( r );
    }

    size_t hard_block_count() const {
        return m_bpr;
    }

    size_t size() const {
        return m_columns * m_rows;
    }

    raw_block_pointer begin_row( size_t idx ) const {
        ASSERT_VALID_RANGE( idx, 0, row_count() )

        return m_data[ idx ].get();
    }

    raw_block_pointer end_row( size_t idx ) const {
        ASSERT_VALID_RANGE( idx, 0, row_count() )
        return m_data[idx].get() + m_data[idx].m_size;
    }

    raw_block_pointer begin_block_row( size_t idx ) const {
        return begin_row( idx );
    }

    raw_block_pointer end_block_row( size_t idx ) const {
        return end_row(idx);
    }

    virtual size_t grow( size_t rows, size_t cols ) {
        this->resize( rows, cols );
        return this->size();
    }

    row_vector & getRow( size_t idx ) {
#ifdef DEBUGGING
        assert( idx < m_rows );
#endif  // DEBUGGING

        return m_data[ idx ];
    }

    void clear( ) {
        m_data.clear();
        m_data.reserve( m_rows );
    }

    virtual ~association_matrix() { }

protected:

/**
 * bpr - blocks_per_row
 */
    inline size_t block_index( size_t row, size_t col, size_t bpr ) const {
        return row * bpr + col / clotho::utility::BitHelper< BlockType >::BITS_PER_BLOCK;
    }

    inline size_t block_row_offset( size_t row, size_t bpr ) const {
        return row * bpr;
    }

    inline size_t block_column_offset( size_t col ) const {
        return col / clotho::utility::BitHelper< BlockType >::BITS_PER_BLOCK;
    }

    void resize( size_t r, size_t c ) {
        m_bpr =  bit_helper_type::padded_block_count( c );
        m_rows = r;
        m_columns = c;

#ifdef DEBUGGING
        std::cerr << "Resizing population to: " << m_rows << " x " << m_columns << " [" << m_bpr << "]" << std::endl;
#endif  // DEBUGGING

        clear();
    }

    data_vector         m_data;    

    size_t m_rows, m_columns;
    size_t m_bpr;
};

}   // namespace genetics
}   // namespace clotho

#undef ASSERT_VALID_RANGE

#endif  // CLOTHO_ROW_GROUPED_ASSOCIATION_MATRIX_HPP_
