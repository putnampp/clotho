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
#ifndef VECTOR_SUBSET_HPP_
#define VECTOR_SUBSET_HPP_

#include "clotho/powerset/config.hpp"

#include <iostream>
#include <algorithm>
#include "clotho/powerset/powerset.hpp"


namespace clotho {
namespace powersets {

template < class Element, class Block = unsigned long, class BlockMap = block_map< Element, Block>, class ElementKeyer = element_key_of< Element > >
class vector_subset {
public:
    typedef vector_subset< Element, Block, BlockMap, ElementKeyer > self_type;
    typedef Element value_type;

    typedef clotho::powersets::powerset< Element, self_type, Block, BlockMap, ElementKeyer > powerset_type;
    typedef std::vector< unsigned int >               data_container_type;

    typedef typename powerset_type::element_index_type  element_index_type;
    typedef std::pair< element_index_type, size_t >     index_type;

    typedef typename data_container_type::iterator          element_index_iterator;
    typedef typename data_container_type::const_iterator    element_index_citerator;

    typedef typename data_container_type::iterator          data_iterator;
    typedef typename data_container_type::const_iterator    data_citerator;

    typedef typename powerset_type::subset_ptr pointer;

    static const element_index_type npos = -1;

    friend class clotho::powersets::powerset< Element, self_type, Block, BlockMap, ElementKeyer >;

    pointer clone() const;

    powerset_type *     getParent() const;

    bool                isSameFamily( pointer other ) const;

    std::pair< typename powerset_type::element_index_type, bool >   addElement( const value_type & elem );
    void                removeElement( const value_type & elem );

    /**
     *  Returns whether an element is a member of this subset
     */
    bool                operator[]( const value_type & elem );

    /**
     *  Returns whether index is in set; equivalent to determining whether
     *  a known element is a member of this subset
     *
     *  Assumes that order_data has already been called
     */
    bool                check_state( element_index_type idx ) const;

    /**
     *  Returns number of blocks allocated for bitset
     */
    size_t              num_blocks() const;

    /**
     *  Returns the number of elements in subset
     */
    size_t              count() const;

    /**
     * Returns the numbers the size of the subset in bits;
     *
     * NOTE: not necessarily equal to the number of total possible elements in the parent set
     */
    size_t              size() const;

    /**
     * Returns the maximum number of elements in parent set
     */
    size_t              max_size() const;

    /**
     *  element_index_iterator for subset
     */
    element_index_iterator begin();
    element_index_iterator end();

    element_index_citerator begin() const;
    element_index_citerator end() const;

    index_type find_first() const;
    index_type find_next( index_type idx ) const;

    /**
     *  Order the data set.
     *
     *  NOTE: This only orders the element indices provided by the parent
     *  powerset.  It does not necessarily reflect an element specific ordering.
     */
    inline void order_data();

    template < class E, class B, class BM, class EK >
    friend bool operator==(const vector_subset< E, B, BM, EK> & lhs, const vector_subset< E, B, BM, EK> & rhs );

    template < class E, class B, class BM, class EK >
    friend std::ostream & operator<<( std::ostream & lhs, const vector_subset< E, B, BM, EK> & rhs );

    virtual ~vector_subset();

    CLOTHO_PROTECTED

    vector_subset(powerset_type * p);
    vector_subset(powerset_type * p, const data_container_type & b);
    vector_subset( const self_type & vs );

    powerset_type *     m_parent;
    data_container_type m_data;
    bool                m_dirty;
};

#define TEMPLATE_HEADER template < class Element, class Block, class BlockMap, class ElementKeyer >
#define SUBSET_SPECIALIZATION vector_subset< Element, Block, BlockMap, ElementKeyer >
#define POWERSET_SPECIALIZATION SUBSET_SPECIALIZATION::powerset_type

TEMPLATE_HEADER
SUBSET_SPECIALIZATION::vector_subset( powerset_type * p ) :
    m_parent(p)
    , m_dirty(false) {
}

TEMPLATE_HEADER
SUBSET_SPECIALIZATION::vector_subset( powerset_type * p, const data_container_type & b ) :
    m_parent(p)
    , m_data( b )
    , m_dirty( true ) {
}

TEMPLATE_HEADER
SUBSET_SPECIALIZATION::vector_subset( const self_type & vs ) :
    m_parent(vs.m_parent)
    , m_data( vs.m_data )
    , m_dirty( vs.m_dirty ) {
}

TEMPLATE_HEADER
SUBSET_SPECIALIZATION::~vector_subset() {}

TEMPLATE_HEADER
typename SUBSET_SPECIALIZATION::powerset_type * SUBSET_SPECIALIZATION::getParent() const {
    return m_parent;
}

TEMPLATE_HEADER
bool SUBSET_SPECIALIZATION::isSameFamily( pointer p ) const {
    return (this->m_parent == p->m_parent);
}

TEMPLATE_HEADER
typename SUBSET_SPECIALIZATION::pointer SUBSET_SPECIALIZATION::clone() const {
    pointer sub = this->m_parent->create_subset( m_data );
    return sub;
}


TEMPLATE_HEADER
bool    SUBSET_SPECIALIZATION::operator[]( const value_type & elem ) {
    order_data();
    typename powerset_type::element_index_type idx = m_parent->find( elem );

    // binary search assumes data is ordered
    return std::binary_search( m_data.begin(), m_data.end(), idx );
}

TEMPLATE_HEADER
std::pair< typename SUBSET_SPECIALIZATION::powerset_type::element_index_type, bool >   SUBSET_SPECIALIZATION::addElement( const value_type & elem ) {
    std::pair< typename powerset_type::element_index_type, bool > res = m_parent->find_or_create( elem );

    if( res.first == npos ) return res;

    m_data.push_back( res.first );

    m_dirty = true;

    return res;
}

TEMPLATE_HEADER
void SUBSET_SPECIALIZATION::removeElement( const value_type & elem ) {
    order_data();

    typename powerset_type::element_index_type idx = m_parent->find( elem );

    if( idx == npos || m_data.size() <= idx ) return;

    // lower_bound assumes data is ordered
    element_index_iterator it = std::lower_bound( m_data.begin(), m_data.end(), idx );

    if( *it == idx ) {
        m_data.erase( it );

        // because data is ordered, erasing an element does not
        // make the list dirty
    }
}

TEMPLATE_HEADER
size_t SUBSET_SPECIALIZATION::count( ) const {
    return m_data.size();
}

TEMPLATE_HEADER
size_t SUBSET_SPECIALIZATION::num_blocks( ) const {
    return m_data.size();
}

TEMPLATE_HEADER
size_t SUBSET_SPECIALIZATION::size() const {
    return m_data.size();
}

TEMPLATE_HEADER
size_t SUBSET_SPECIALIZATION::max_size() const {
    return m_parent->variable_allocated_size();
}

TEMPLATE_HEADER
bool SUBSET_SPECIALIZATION::check_state( element_index_type idx ) const {
    assert( !m_dirty );
    return std::binary_search( m_data.begin(), m_data.end(), idx );
}

TEMPLATE_HEADER
typename SUBSET_SPECIALIZATION::element_index_iterator SUBSET_SPECIALIZATION::begin() {
    return m_data.begin();
}

TEMPLATE_HEADER
typename SUBSET_SPECIALIZATION::element_index_iterator SUBSET_SPECIALIZATION::end() {
    return m_data.end();
}

TEMPLATE_HEADER
typename SUBSET_SPECIALIZATION::element_index_citerator SUBSET_SPECIALIZATION::begin() const {
    return m_data.begin();
}

TEMPLATE_HEADER
typename SUBSET_SPECIALIZATION::element_index_citerator SUBSET_SPECIALIZATION::end() const {
    return m_data.end();
}

TEMPLATE_HEADER
typename SUBSET_SPECIALIZATION::index_type SUBSET_SPECIALIZATION::find_first() const {
    if( m_data.empty() ) {
        return std::make_pair( (element_index_type) npos, (size_t) npos );
    }
    return std::make_pair( m_data[0], (size_t) 0);
}

TEMPLATE_HEADER
typename SUBSET_SPECIALIZATION::index_type SUBSET_SPECIALIZATION::find_next( index_type idx ) const {
    if( idx.second == (size_t)npos)
        return std::make_pair( (element_index_type) npos, (size_t) npos );

    size_t tmp = idx.second + 1;
    if( tmp >= m_data.size() ) {
        return std::make_pair( (element_index_type) npos, (size_t) npos );
    } else {
        return std::make_pair( m_data[tmp], tmp );
    }
}

TEMPLATE_HEADER
void SUBSET_SPECIALIZATION::order_data() {
    if( m_dirty ) {
        std::sort( m_data.begin(), m_data.end() );
        m_dirty = false;
    }
}

/**
 * Equality of vector_subsets.
 */
TEMPLATE_HEADER
inline bool operator==( const SUBSET_SPECIALIZATION & lhs, const SUBSET_SPECIALIZATION & rhs ) {
    if( lhs.m_parent != rhs.m_parent || lhs.m_data.size() != rhs.m_data.size() )
        return false;

    // if one of the  sets is dirty we have to assume
    // that elements are out of order;
    // This is expensive, but easily avoided if order_data() is called
    // before equality is checked.
    // Therefore, will just assert if either is dirty, in the hopes
    // that the developer will simply add the order_data call before
    // their check statement.
    //
    assert( !lhs.m_dirty && !rhs.m_dirty );

    return (lhs.m_data == rhs.m_data);
}

TEMPLATE_HEADER
std::ostream & operator<<( std::ostream & lhs, const SUBSET_SPECIALIZATION & rhs ) {
    lhs << "{<";

    if( !rhs.m_data.empty() ) {
        typename SUBSET_SPECIALIZATION::element_index_citerator it = rhs.begin();
        lhs << *it;
        while( ++it != rhs.end() ) {
            lhs << "," << *it;
        }
    }
    lhs << ">; " << rhs.m_data.size() << "}";
    return lhs;
}

#undef SUBSET_SPECIALIZATION
#undef TEMPLATE_HEADER
#undef POWERSET_SPECIALIZATION

}   // namespace powersets
}   // namespace clotho

#include "clotho/powerset/vector_subset_iterator.hpp"

#endif  // VECTOR_SUBSET_HPP_
