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
#ifndef CLOTHO_ASSOCIATION_MATRIX_HPP_
#define CLOTHO_ASSOCIATION_MATRIX_HPP_

#include "clotho/data_spaces/growable2D.hpp"
#include "clotho/utility/bit_helper.hpp"
#include <cstring>
//#include <iostream>

namespace clotho {
namespace genetics {

#ifdef DEBUGGING

#define ASSERT_VALID_RANGE( x, min, max) assert( min <= x && x < max)

#else   // DEBUGGING

#define ASSERT_VALID_RANGE( x, min, max)

#endif  // DEBUGGING

template < class BlockType = unsigned long long >
class association_matrix : public growable2D {
public:
    typedef BlockType   block_type;

    typedef clotho::utility::BitHelper< block_type >    bit_helper_type;

    class BlockIterator {
    public:
        BlockIterator( block_type * start, size_t row_bound, size_t column_bound, size_t step = 1 ) : 
            m_start( start )
            , m_ptr( start )
            , m_row_bound( row_bound )
            , m_column_bound( column_bound )
            , m_step( step )
            , m_i(0)
            , m_j(0)
        {}

        BlockIterator( const BlockIterator & other ) : 
            m_start( other.m_start )
            , m_ptr( other.m_ptr )
            , m_row_bound( other.m_row_bound )
            , m_column_bound( other.m_column_bound )
            , m_step( other.m_step )
            , m_i( other.m_i )
            , m_j( other.m_j )
        {}

        bool hasNext() const {
            return ( m_i < m_row_bound );
        }

        block_type next() {
            assert( hasNext() );

            block_type v = *m_ptr;

            advance_pointer();
            
            return v;
        }

        std::pair< block_type, block_type > next_pair() {
            assert( hasNext() );

            block_type v = *m_ptr++;
            block_type w = *m_ptr;

            advance_pointer();
            
            return std::make_pair(v, w);
        }

        void write_next( block_type v ) {
            assert( hasNext() );

            *(m_ptr) = v;
            advance_pointer();
        }

        size_t max_rows() const {
            return m_row_bound;
        }

        size_t max_columns() const {
            return m_column_bound;
        }

        virtual ~BlockIterator() {}
    protected:

        void advance_pointer() {
            m_j += m_step;
            if( m_j >= m_column_bound ) {
                m_j = 0;
                ++m_i;
            }

            m_ptr = m_start + m_i * m_column_bound + m_j;
        }

        block_type      * m_start, *m_ptr;
        size_t          m_row_bound, m_column_bound, m_step;
        size_t          m_i, m_j;
    };

    typedef BlockIterator       block_iterator;
    typedef BlockIterator       row_iterator;
    typedef BlockIterator       column_iterator;
    typedef BlockIterator       row_pair_iterator;


    association_matrix( size_t rows = 1, size_t columns = 1 ) :
        m_data( NULL )
        , m_rows(0)
        , m_columns(0)
        , m_allocated_size(0)
    {
        this->resize( rows, columns );
    }

    bool operator()( size_t row, size_t col ) {
        ASSERT_VALID_RANGE( row, 0, m_rows );
        ASSERT_VALID_RANGE( col, 0, m_columns );

        return ((m_data[ block_index( row, col ) ] & bit_helper_type::bit_offset( col )) != 0);
    }

    void flip( size_t row, size_t col ) {
        ASSERT_VALID_RANGE( row, 0, m_rows );
        ASSERT_VALID_RANGE( col, 0, m_columns );

        m_data[ block_index( row, col ) ] ^= bit_helper_type::bit_offset( col );
    }

    size_t row_count() const {
        return m_rows;
    }

    size_t column_count() const {
        return m_columns;
    }

    size_t block_row_count() const {
        return scale_columns( m_columns );
    }

    size_t block_column_count() const {
        return m_rows;
    }

    size_t block_count() const {
        return m_rows * scale_columns( m_columns );
    }

    block_iterator raw_iterator() const {
        return block_iterator( m_data, block_row_count(), block_column_count() );
    }

    row_iterator   getRowAt( size_t idx ) const {
        ASSERT_VALID_RANGE( idx, 0, block_column_count() );

        return row_iterator( m_data + idx, block_row_count(), block_column_count(), block_column_count() );
    }

    row_pair_iterator getRowPairAt( size_t idx ) const {
        ASSERT_VALID_RANGE( idx, 0, block_column_count() );
        ASSERT_VALID_RANGE( idx + 1, 0, block_column_count() );
        return row_iterator( m_data + idx, block_row_count(), block_column_count(), block_column_count() );
    }

    column_iterator getColumnAt( size_t idx ) const {
        ASSERT_VALID_RANGE( idx, 0, block_row_count() );

        return column_iterator( m_data + idx * block_column_count(), 1, block_column_count() );
    }

    size_t size() const {
        return m_columns * m_rows;
    }

    size_t allocated_size() const {
        return m_allocated_size;
    }

    virtual size_t grow( size_t rows, size_t cols ) {
        this->resize( rows, cols );
        return this->size();
    }

    void clear() {
        memset( m_data, 0, sizeof( block_type ) * m_allocated_size);
    }

    virtual ~association_matrix() {
        if( m_data != NULL ) {
            delete [] m_data;
        }
    }

protected:
    
    inline size_t block_index( size_t row, size_t col ) const {
        return ((col / bit_helper_type::BITS_PER_BLOCK) * m_rows) + row;
    }

    size_t scale_columns( size_t cols ) const {
        return (cols / bit_helper_type::BITS_PER_BLOCK + 1 );
    }

    void resize( size_t rows, size_t columns ) {
        size_t block_rows = scale_columns( columns );
        size_t new_size = rows * block_rows;

        if( new_size > m_allocated_size ) {
            if( m_data != NULL ) {
                delete [] m_data;
            }

            m_data = new block_type[ new_size ];
            m_allocated_size = new_size;
        }
        
        m_rows = rows;
        m_columns = columns;
    }

    block_type  * m_data;

    size_t m_rows, m_columns, m_allocated_size;
};

}   // namespace genetics
}   // namespace clotho

#undef ASSERT_VALID_RANGE

#endif  // CLOTHO_ASSOCIATION_MATRIX_HPP_
