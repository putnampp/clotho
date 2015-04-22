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
#ifndef CLOTHO_VECTOR_SUBSET_ITERATOR_HELPER_HPP_
#define CLOTHO_VECTOR_SUBSET_ITERATOR_HELPER_HPP_

#include "clotho/utility/iterator_helper.hpp"
#include "clotho/powerset/vector_subset.hpp"

namespace clotho {
namespace utility {

template < class E, class B, class BK, class EK >
class element_iterator< clotho::powersets::vector_subset< E, B, BK, EK > > {
public:
    typedef clotho::powersets::vector_subset< E, B, BK, EK >    container_type;
    typedef element_iterator< container_type >                  self_type;
    typedef E                                                   value_type;

    typedef typename container_type::data_container_type::size_type size_type;

    typedef typename container_type::powerset_type::variable_iterator  iterator;
    typedef typename container_type::element_index_iterator           index_iterator;

    static const size_t npos = -1;

    element_iterator( container_type & con ) :
        m_container( &con )
        , m_idx( con.begin() )
        , m_var(con.getParent()->variable_begin())
        , m_offset( npos ) {
        m_offset = ((m_idx != con.end() ) ? *m_idx : m_offset);
    }

    element_iterator( container_type & con, index_iterator idx ) :
        m_container( &con )
        , m_idx( idx )
        , m_var( con.getParent()->variable_begin() )
        , m_offset( npos ) {
        m_offset = ((m_idx != con.end() ) ? *m_idx : m_offset);
    }

    element_iterator( const self_type & other ) :
        m_container( other.m_container )
        , m_idx( other.m_idx )
        , m_var( other.m_var )
        , m_offset( other.m_offset ) {
    }

    self_type & operator++() {
        if( m_idx != m_container->end() ) {
            m_offset = (( ++m_idx == m_container->end() ) ? npos : *m_idx);
        }
        return *this;
    }

    self_type operator++( int ) {
        self_type tmp(*this);
        operator++();
        return tmp;
    }

    value_type &  operator*() {
        assert( m_offset != npos );

        return (*(m_var + m_offset));
    }

    bool operator==( const self_type & rhs ) const {
        return (m_container == rhs.m_container && m_idx == rhs.m_idx);
    }

    bool operator!=( const self_type & rhs ) const {
        return (m_container != rhs.m_container || m_idx != rhs.m_idx);
    }

    virtual ~element_iterator() {}
protected:
    container_type  *   m_container;
    index_iterator      m_idx;
    iterator            m_var;

    size_t              m_offset;
};

template < class E, class B, class BK, class EK >
struct iterator_helper< clotho::powersets::vector_subset< E, B, BK, EK > > {
    typedef clotho::powersets::vector_subset< E, B, BK, EK > subset_type;
    typedef element_iterator< subset_type > type;

    static type make_first( subset_type & s ) {
        return type( s );
    }

    static type make_last( subset_type & s ) {
        return type( s, s.end() );
    }

    static std::pair< type, type > make_range( subset_type & s ) {
        return std::make_pair( type(s), type(s, s.end()) );
    }
};

}   // namespace utility {
}   // namespace clotho {

#endif  // CLOTHO_VECTOR_SUBSET_ITERATOR_HELPER_HPP_
