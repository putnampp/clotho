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
#include <deque>

#include "clotho/data_spaces/growable1D.hpp"

#include <iostream>

namespace clotho {
namespace genetics {

template < class PositionType = double >
class base_allele_vectorized : public growable1D {
public:
    typedef base_allele_vectorized< PositionType >              self_type;
    typedef PositionType                                        position_type;

    typedef std::vector< position_type >                        position_vector_type;
    typedef typename position_vector_type::iterator             position_iterator;
    typedef typename position_vector_type::const_iterator       const_position_iterator;

    typedef std::deque< size_t >                                free_vector_type;
    typedef typename free_vector_type::iterator                 free_iterator;
    typedef typename free_vector_type::const_iterator           const_free_iterator;

    base_allele_vectorized( size_t a = 0 ) {
        this->grow( a );
    }

    position_type getPositionAt( size_t index ) const {
        return m_positions[ index ];
    }

    void setPositionAt( size_t index, position_type pos ) {
        if( index >= m_positions.size() ) {
            resize( index + 1 );
        }

        m_positions[ index ] = pos;
    }

    position_iterator position_begin() {
        return m_positions.begin();
    }

    position_iterator position_end() {
        return m_positions.end();
    }

    const_position_iterator position_begin() const {
        return m_positions.begin();
    }

    const_position_iterator position_end() const {
        return m_positions.begin();
    }

    void inherit( self_type & parent ) {
        std::copy( parent.position_begin(), parent.position_end(), this->position_begin() );
    }

    size_t next_free() {
        if( m_free.empty() ) {
            return -1;
        }

        size_t res = m_free.front();
        m_free.pop_front();
        return res;
    }

    virtual size_t grow( size_t s ) {
        resize( s );

        return size();
    }

    size_t size() const {
        return m_positions.size();
    }

    template < class Iterator >
    void updateFreeSpace( Iterator start, Iterator end ) {
        while( start != end ) {
            m_free.push_back( *start );
            ++start;
        }
    }

    size_t      free_size() const {
        return m_free.size();
    }

    free_iterator free_begin() {
        return m_free.begin();
    }

    free_iterator free_end() {
        return m_free.end();
    }
    
    const_free_iterator free_begin() const {
        return m_free.begin();
    }

    const_free_iterator free_end() const {
        return m_free.end();
    }

    void push_back( self_type & other, size_t idx ) {
        size_t e = this->m_positions.size();

        this->resize( e + 1 );

        this->setPositionAt( e, other.getPositionAt( idx ) );
    }

    virtual ~base_allele_vectorized() {}

protected:

    virtual void resize( size_t s ) {
        size_t n = m_positions.size();

        m_positions.resize( s );

        size_t m = m_positions.size();

        while( n < m ) {
            m_free.push_back( n );
            ++n;
        }
    }

    position_vector_type            m_positions;

    free_vector_type                m_free;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_BASE_ALLELE_HPP_
