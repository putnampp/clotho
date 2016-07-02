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
#ifndef CLOTHO_ASSOCIATION_VECTOR_HPP_
#define CLOTHO_ASSOCIATION_VECTOR_HPP_

#include <memory>
#include <algorithm>
#include "clotho/utility/bit_helper.hpp"

#ifdef DEBUGGING
#define ASSERT_LTE(x, y) assert( x <= y );
#else
#define ASSERT_LTE(x, y)
#endif  // DEBUGGING

namespace clotho {
namespace genetics {

template < class BlockType >
struct association_vector {
    typedef association_vector< BlockType > self_type;

    typedef BlockType                       block_type;
    typedef clotho::utility::BitHelper< block_type > bit_helper_type;
    typedef std::shared_ptr< block_type >   pointer;

    typedef block_type *                    raw_pointer;
    typedef const block_type *              const_raw_pointer;

    size_t  m_size;
    pointer m_data;
    bool    m_readonly;

    association_vector( size_t N ) :
        m_size( N )
        , m_data()
        , m_readonly(false)
    {
        assert( N > 0);
        m_data = pointer( new block_type[ N ], std::default_delete< block_type[] >() );
        raw_pointer p = m_data.get();
        for( size_t i = 0; i < N; ++i )
            p[ i ] = bit_helper_type::ALL_UNSET;
    }

    association_vector( const self_type & other ) :
        m_size(other.m_size)
        , m_data( other.m_data )
        , m_readonly( other.m_readonly )
    {    }

    association_vector( size_t N, raw_pointer ptr, bool readonly = false ) :
        m_size( N )
        , m_data( ptr, std::default_delete< block_type[] >() )
        , m_readonly( readonly )
    {}

    raw_pointer get() {
        return m_data.get();
    }

    const_raw_pointer get() const {
        return m_data.get();
    }

    void copy( const self_type & other ) {
        ASSERT_LTE( other.m_size, m_size )
        std::copy( other.get(), other.get() + other.m_size, m_data.get() );
    }

    bool operator==( const self_type & rhs ) {
#ifdef DEBUGGING
        if( this->m_data == rhs.m_data) {
            assert( this->m_readonly == rhs.m_readonly && this->m_size == rhs.m_size );
        }
#endif  //DEBUGGING
        return this->m_data == rhs.m_data;
    }

    void finalize() {
        m_readonly = true;
    }

    virtual ~association_vector() {}
};

}   // namespace genetics
}   // namespace clotho

#ifdef DEBUGGING
#undef ASSERT_LTE
#endif  // DEBUGGING

#endif  // CLOTHO_ASSOCIATION_VECTOR_HPP_
