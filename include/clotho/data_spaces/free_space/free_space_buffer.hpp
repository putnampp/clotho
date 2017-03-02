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
#ifndef CLOTHO_FREE_SPACE_BUFFER_HPP_
#define CLOTHO_FREE_SPACE_BUFFER_HPP_

#include "clotho/utility/bit_helper.hpp"

namespace clotho {
namespace genetics {

template < class BlockType >
class free_space_buffer {
public:
    typedef free_space_buffer< BlockType > self_type;
    typedef BlockType block_type;
    typedef clotho::utility::BitHelper< BlockType > bit_helper_type;

    free_space_buffer( ) :
        m_var_buffer( NULL )
        , m_fixed_buffer( NULL )
        , m_temp_buffer( NULL )
        , m_buffer_size(0)
        , m_buffer_alloc(0)
    {}

    void updateBuffers( block_type * buf, const unsigned int len ) {
        assert( len <= m_buffer_size);

        updateVariableBuffer( buf, len );
        updateFixedBuffer( buf, len );
    }

    template < class Iterator >
    void updateBuffers( Iterator first, Iterator last ) {
        memset( m_temp_buffer, 0, sizeof(block_type) * m_buffer_size);
        while( first != last ) {
            unsigned int idx = (*first / bit_helper_type::BITS_PER_BLOCK);
            block_type b = ((block_type) 1 << (*first % bit_helper_type::BITS_PER_BLOCK ));

            m_temp_buffer[ idx ] |= b;
            ++first;
        }

        updateBuffers( m_temp_buffer, m_buffer_size );
    }

    block_type * getVariableBuffer() {
        return m_var_buffer;
    }

    block_type * getFixedBuffer() {
        return m_fixed_buffer;
    }

    unsigned int size() const {
        return m_buffer_size;
    }

    void reset( unsigned int s ) {
        if( m_buffer_alloc < s ) {
            if( m_var_buffer != NULL ) {
                delete [] m_var_buffer;
            }

            if( m_fixed_buffer != NULL ) {
                delete [] m_fixed_buffer;
            }

            if( m_temp_buffer != NULL ) {
                delete [] m_temp_buffer;
            }

            m_var_buffer = new block_type[ s ];
            m_fixed_buffer = new block_type[ s ];
            m_temp_buffer = new block_type[ s ];

            m_buffer_alloc = s;
        }

        m_buffer_size = s;
        memset( m_var_buffer, 0, m_buffer_alloc * sizeof(block_type) );
        memset( m_temp_buffer, 0, m_buffer_alloc * sizeof(block_type) );
        memset( m_fixed_buffer, 255, m_buffer_alloc * sizeof(block_type) );
    }

    virtual ~free_space_buffer() {
        if( m_var_buffer != NULL ) {
            delete [] m_var_buffer;
        }

        if( m_fixed_buffer != NULL ) {
            delete [] m_fixed_buffer;
        }

        if( m_temp_buffer != NULL ) {
            delete [] m_temp_buffer;
        }
    }

protected:

    void updateVariableBuffer( block_type * buf, const unsigned int len ) {
        unsigned int i = 0;
        while( i < len ) {
            m_var_buffer[i] |= buf[i];
            ++i;
        }

    }

    void updateFixedBuffer( block_type * buf, const unsigned int len ) {
        unsigned int i = 0;
        while( i < len ) {
            m_fixed_buffer[ i ] &= buf[i];
            ++i;
        }
    }

    block_type * m_var_buffer;
    block_type * m_fixed_buffer;
    block_type * m_temp_buffer;

    unsigned int m_buffer_size, m_buffer_alloc;
};

}   // namespace genetics
}   // namespace clotho
#endif  // CLOTHO_FREE_SPACE_BUFFER_HPP_
