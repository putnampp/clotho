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
#ifndef CLOTHO_FREE_SPACE_ACCUMULATOR_HPP_
#define CLOTHO_FREE_SPACE_ACCUMULATOR_HPP_

#include "clotho/data_spaces/free_space/free_space_buffer.hpp"
#include "clotho/data_spaces/free_space/free_space_details.hpp"

namespace clotho {
namespace genetics {

template < class BlockType, class SizeType >
class free_space_accumulator : public free_space_details< SizeType >, public free_space_buffer< BlockType > {
public:
    typedef free_space_accumulator< BlockType, SizeType >   self_type;
    typedef free_space_details< SizeType >                  base_type;
    typedef free_space_buffer< BlockType >                  buffer_type;

    free_space_accumulator( ) {}

    void resetBuffers( buffer_type * other ) {
        this->reset( other->size() );
        updateBuffers( other );
    }

    void updateBuffers( buffer_type * buf ) {
        this->updateVariableBuffer( buf->getVariableBuffer(), buf->size());
        this->updateFixedBuffer( buf->getFixedBuffer(), buf->size() );
    }

    void finalize( unsigned int max_alleles ) {
        this->resize( max_alleles );
        this->analyze_free_indices( this->getFixedBuffer(), this->getVariableBuffer(), this->size(), max_alleles );
    }

    virtual ~free_space_accumulator() { }

protected:
};


template < class BlockType, class SizeType >
class free_space_accumulator_mt : public free_space_details< SizeType > {
public:
    typedef free_space_accumulator_mt< BlockType, SizeType >   self_type;
    typedef free_space_details< SizeType >                  base_type;

    typedef BlockType                                       block_type;

    typedef clotho::utility::BitHelper< block_type >        bit_helper_type;

    struct thread_buffer {
        free_space_accumulator_mt * m_parent;
        block_type * m_thread_var_buffer;
        block_type * m_thread_fix_buffer;

        thread_buffer() :
            m_parent( NULL )
            , m_thread_var_buffer( NULL )
            , m_thread_fix_buffer( NULL )
        {}

        thread_buffer( free_space_accumulator_mt * par, block_type * var_buf, block_type * fix_buf ) :
            m_parent( par )
            , m_thread_var_buffer( var_buf )
            , m_thread_fix_buffer( fix_buf )
        {}

        thread_buffer( const thread_buffer & other ) :
            m_parent( other.m_parent )
            , m_thread_var_buffer( other.m_thread_var_buffer )
            , m_thread_fix_buffer( other.m_thread_fix_buffer )
        {}
            
        void update( block_type * buf, unsigned int len ) {
            assert( m_parent != NULL );

            m_parent->updateVariableBuffer( m_thread_var_buffer, buf, len );
            m_parent->updateFixedBuffer( m_thread_fix_buffer, buf, len );
        }

        void reset( const thread_buffer & other ) {
            m_parent = other.m_parent;
            m_thread_var_buffer = other.m_thread_var_buffer;
            m_thread_fix_buffer = other.m_thread_fix_buffer;
        }

        virtual ~thread_buffer() {}
    };

    typedef thread_buffer                                   buffer_type;

    free_space_accumulator_mt( unsigned int tc ) :
        m_var_buffer( NULL )
        , m_fix_buffer( NULL )
        , m_thread_count( tc )
        , m_blocks_per_row(0)
        , m_buffer_alloc(0)
    {}

    void resetBuffers( unsigned int block_count ) {
        if( block_count > m_blocks_per_row ) {
            if( m_var_buffer != NULL ) {
                delete [] m_var_buffer;
            }

            if( m_fix_buffer != NULL ) {
                delete [] m_fix_buffer;
            }

            m_buffer_alloc = m_thread_count * block_count;

            m_var_buffer = new block_type[ m_buffer_alloc ];
            m_fix_buffer = new block_type[ m_buffer_alloc ];
        }

        m_blocks_per_row = block_count;

        memset( m_var_buffer, 0, m_buffer_alloc * sizeof(block_type));
        memset( m_fix_buffer, 255, m_buffer_alloc * sizeof(block_type));
    }

    buffer_type getThreadBuffer( unsigned int i ) {
        assert( i < m_thread_count );

        i *= m_blocks_per_row;

        return buffer_type( this, m_var_buffer + i, m_fix_buffer + i );
    }

    void finalize( unsigned int max_alleles ) {
        finalize( max_alleles, m_thread_count );
    }

    void finalize( unsigned int max_alleles, unsigned int T ) {
        this->resize( max_alleles );
        this->reduce( T );

        this->analyze_free_indices( m_fix_buffer, m_var_buffer, m_blocks_per_row, max_alleles );
    }

    void updateVariableBuffer( block_type * dest, block_type * source, unsigned int N ) {
        assert( N <= m_blocks_per_row);
        while( N-- ) {
            *dest++ |= *source++;
        }
    }

    void updateFixedBuffer( block_type * dest, block_type * source, unsigned int N ) {
        assert( N <= m_blocks_per_row);
        unsigned int tail = m_blocks_per_row - N;
        while( N-- ) {
            *dest++ &= *source++;
        }

        // clear any trailing blocks
        while( tail-- ) {
            *dest++ = bit_helper_type::ALL_UNSET;
        }
    }

    virtual ~free_space_accumulator_mt() {
        if( m_var_buffer != NULL ) {
            delete [] m_var_buffer;
        }

        if( m_fix_buffer != NULL ) {
            delete [] m_fix_buffer;
        }
    }

protected:

    void reduce( unsigned int T ) {
        assert( T <= m_thread_count );

        for( unsigned int i = 1, offset = m_blocks_per_row; i < T; ++i, offset += m_blocks_per_row ) {
            updateVariableBuffer( m_var_buffer, m_var_buffer + offset, m_blocks_per_row );
            updateFixedBuffer( m_fix_buffer, m_fix_buffer + offset, m_blocks_per_row );
        }
    }

    block_type * m_var_buffer;
    block_type * m_fix_buffer;
    unsigned int m_thread_count, m_blocks_per_row, m_buffer_alloc;
};

}   // namespace genetics
}   // namespace clotho

namespace clotho {
namespace utility {

template < class BlockType, class SizeType >
struct state_getter< clotho::genetics::free_space_accumulator< BlockType, SizeType > > {
    typedef clotho::genetics::free_space_accumulator< BlockType, SizeType > object_type;

    void operator()( boost::property_tree::ptree & s, object_type & obj ) {
        boost::property_tree::ptree fr;
        clotho::utility::add_value_array( fr, obj.free_begin(), obj.free_end() );

        s.put_child( "free", fr );
        s.add( "free_size", obj.free_size() );
        s.add( "variable_size", obj.variable_count() );
    }
};

template < class BlockType, class SizeType >
struct state_getter< clotho::genetics::free_space_accumulator_mt< BlockType, SizeType > > {
    typedef clotho::genetics::free_space_accumulator_mt< BlockType, SizeType > object_type;

    void operator()( boost::property_tree::ptree & s, object_type & obj ) {
        boost::property_tree::ptree fr;
        clotho::utility::add_value_array( fr, obj.free_begin(), obj.free_end() );

        s.put_child( "free", fr );
        s.add( "free_size", obj.free_size() );
        s.add( "variable_size", obj.variable_count() );
    }
};
}   // namespace utility
}   // namespace clotho
#endif  // CLOTHO_FREE_SPACE_ACCUMULATOR_HPP_
