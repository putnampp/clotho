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
#ifndef CLOTHO_POPULATION_MATRIX_HPP_
#define CLOTHO_POPULATION_MATRIX_HPP_

#include <utility>

#define _DEBUGGING_

#if _DEBUGGING_

#include <cassert>
#define INDEX_CHECK( idx ) assert( 0 < idx && idx < m_height );
#define ALLOCATION_CHECK( ptr ) assert( ptr != NULL );

#else

#define INDEX_CHECK( idx )
#define ALLOCATION_CHECK( ptr )

#endif  // _DEBUGGING_

namespace clotho {
namespace genetics {

template < class BinaryType = unsigned long long, class SizeType = size_t >
class PopulationMatrix {
public:
    typedef BinaryType  bit_block_type;
    typedef SizeType    size_type;

    typedef bit_block_type * row_pointer;

    static const unsigned int BITS_PER_BLOCK = sizeof( bit_block_type ) * 8; 

    PopulationMatrix( size_type height, size_type bit_width ) :
        m_data(NULL)
        , m_soft_max(NULL)
        , m_height(0)
        , m_width(0)
        , m_size(0)
        , m_soft_size(0)
    {
        resize( height, bit_width );
    }

    row_pointer getRowStart( size_type index ) {
        INDEX_CHECK( index )
        return m_data + (index) * m_width;
    }

    row_pointer getRowEnd( size_type index ) {
        INDEX_CHECK( index )
        return m_data + (index + 1) * m_width;
    }

    std::pair< row_pointer, row_pointer > getRowRange( size_type index ) {
        INDEX_CHECK( index )

        row_pointer s = m_data + ( index * m_width );
        row_pointer e = s + m_width;

        return std::make_pair(s, e);
    }

    size_type getSoftMax( size_t index ) const {
        INDEX_CHECK( index )

        return m_soft_max[ index ];
    }

    size_type getHardMax( ) const {
        return getBitWidth();
    }

    size_type getHeight() const {
        return m_height;
    }

    size_type getBlockWidth() const {
        return m_width;
    }

    size_type getBitWidth() const {
        return m_width * BITS_PER_BLOCK;
    }

    virtual ~PopulationMatrix() {
        if( m_data != NULL ) {
            delete [] m_data;
        }

        if( m_soft_max != NULL ) {
            delete [] m_soft_max;
        }
    }

protected:

    size_type scale_bits_to_blocks( size_type bits ) {
        return (bits / BITS_PER_BLOCK ) + 1;
    }

    void resize( size_type height, size_type bit_width ) {
        size_type _width = scale_bits_to_blocks( bit_width );

        size_type new_size = height * _width;

        if( new_size > m_size ) {

            if( m_data != NULL ) {
                delete [] m_data;
            }

            m_data = new bit_block_type[ new_size ];

            ALLOCATION_CHECK( m_data )

            m_size = new_size;
        }

        if( height > m_soft_size ) {

            if( m_soft_size != NULL ) {
                delete [] m_soft_size;
            }

            m_soft_size = new size_type[ height ];

            m_soft_size = height;
        }

        m_height = height;
        m_width = _width;
    }

    bit_block_type  * m_data;
    soft_size_type  * m_soft_max;
    size_type        m_height, m_width, m_size, m_soft_size;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_POPULATION_MATRIX_HPP_
