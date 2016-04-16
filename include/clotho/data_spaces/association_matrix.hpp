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

//    class BlockIterator {
//    public:
//        typedef block_type  value_type;
//
//        BlockIterator( block_type * start, block_type * end, size_t offset ) : 
//            m_start( start )
//            , m_end( end )
//            , m_offset( offset )
//        {}
//
//        BlockIterator( const BlockIterator & other ) : 
//            m_start( other.m_start )
//            , m_end( other.m_end )
//            , m_offset( other.m_offset )
//        {}
//
//        bool hasNext() const {
//            return (m_start + offset) < m_end;
//        }
//
//        value_type next() {
//            m_start += offset;
//            return *m_start;
//        }
//
//        virtual ~BlockIterator() {}
//    protected:
//        block_type *    m_start, * m_end;
//        
//    };
//
//    typedef BlockIterator       block_iterator;
//

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

    virtual ~association_matrix() {
        if( m_data != NULL ) {
            delete [] m_data;
        }
    }

protected:
    
    inline size_t block_index( size_t row, size_t col ) {
        return ((col / bit_helper_type::BITS_PER_BLOCK) * m_rows) + row;
    }

    size_t scale_columns( size_t cols ) {
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
