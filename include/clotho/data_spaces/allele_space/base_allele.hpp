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
#ifndef CLOTHO_BASE_ALLELE_HPP_
#define CLOTHO_BASE_ALLELE_HPP_

#include <vector>
//#include <deque>

#include "clotho/data_spaces/growable1D.hpp"

#include <iostream>

#include "clotho/utility/state_object.hpp"

#include <boost/dynamic_bitset.hpp>

namespace clotho {
namespace genetics {

template < class PositionType = double >
class base_allele_vectorized : public growable1D {
public:
    typedef base_allele_vectorized< PositionType >              self_type;
    typedef PositionType                                        position_type;

    typedef position_type *                                     position_vector_type;
    typedef position_type *                                     position_iterator;

    typedef boost::dynamic_bitset< unsigned int >               free_vector_type;

    base_allele_vectorized( size_t a = 0 ) :
        m_positions( NULL )
        , m_pos_size(0)
        , m_alloc_size(0)
    {
        this->grow( a );
    }

    position_vector_type getPositions() {
        return m_positions;
    }

    position_type getPositionAt( size_t index ) const {
        return m_positions[ index ];
    }

    void setPositionAt( size_t index, position_type pos ) {
        if( index >= size() ) {
            resize( index + 1 );
        }

        m_positions[ index ] = pos;
    }

    position_iterator position_begin() {
        return m_positions;
    }

    position_iterator position_end() {
        return m_positions + m_pos_size;
    }

    void inherit( self_type & parent ) {
        assert( parent.size() <= size() );

        memcpy( m_positions, parent.m_positions, parent.size() * sizeof(position_type) );
    }

    size_t next_free() {
        size_t res = m_free.find_first();
        if( res != free_vector_type::npos ) {
            m_free.reset( res );
        }
        return res;
    }

    virtual size_t grow( size_t s ) {
        resize( s );

        return size();
    }

    size_t size() const {
        return m_pos_size;
    }

    void trimFreeSpace( size_t parent_size ) {
        while(parent_size < size() ) {
            m_free.set( parent_size++ );
        }       
    }

    template < class Iterator >
    void updateFreeSpace( Iterator start, Iterator end ) {
        while( start != end ) {
            m_free.set( *start++);
        }
    }

    void push_back( self_type & other, size_t idx ) {
        size_t e = this->m_positions.size();

        this->resize( e + 1 );

        this->setPositionAt( e, other.getPositionAt( idx ) );
    }

    virtual ~base_allele_vectorized() {}

protected:

    virtual void resize( size_t s ) {
        if ( s > m_alloc_size ) {
            if( m_positions != NULL ) {
                delete [] m_positions;
            }

            m_positions = new position_type[ s ];

            m_alloc_size = s;
        }

        m_pos_size = s;
        m_free.resize( s );
        m_free.reset();
    }

    position_vector_type            m_positions;
    free_vector_type                m_free;

    size_t  m_pos_size, m_alloc_size;
};

}   // namespace genetics
}   // namespace clotho

namespace clotho {
namespace utility {

template < class PositionType >
struct state_getter< clotho::genetics::base_allele_vectorized< PositionType > > {
    typedef clotho::genetics::base_allele_vectorized< PositionType > object_type;

    void operator()( boost::property_tree::ptree & s, object_type & obj ) {

        size_t all_count = obj.allele_count();
        size_t i = 0;
        while( i < all_count ) {
            boost::property_tree::ptree all;
            all.put( "position", obj.getPositionAt(i) );
            s.push_back( std::make_pair("", all ) );
            ++i;
        }

    }
};

}   // namespace utility
}   // namespace clotho

#endif  // CLOTHO_BASE_ALLELE_HPP_
