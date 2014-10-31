#ifndef VARIABLE_SUBSET_HPP_
#define VARIABLE_SUBSET_HPP_

#include "clotho/powerset/config.hpp"

#include <iostream>
#include "clotho/powerset/powerset.hpp"

namespace clotho {
namespace powersets {

template < class Element, class Block = unsigned long, class BlockMap = block_map< Element, Block>, class ElementKeyer = element_key_of< Element > >
class variable_subset {
public:
    typedef variable_subset< Element, Block, BlockMap, ElementKeyer > self_type;
    typedef Element value_type;

    typedef clotho::powersets::powerset< Element, variable_subset< Element, Block, BlockMap, ElementKeyer >, Block, BlockMap, ElementKeyer > powerset_type;
    typedef boost::dynamic_bitset< unsigned long > bitset_type;
    typedef size_t index_type;

    typedef std::vector< bitset_type::block_type, bitset_type::allocator_type >::iterator block_iterator;
    typedef std::vector< bitset_type::block_type, bitset_type::allocator_type >::const_iterator cblock_iterator;

    typedef typename powerset_type::subset_ptr pointer;

    friend class clotho::powersets::powerset< Element, variable_subset< Element, Block, BlockMap, ElementKeyer >, Block, BlockMap, ElementKeyer >;

    pointer clone() const;

//    void    copy();
//    void                release();
//
//    unsigned int        ref_count() const;

    powerset_type *     getParent() const;

    bool                isSameFamily( pointer other ) const;

    void                addElement( const value_type & elem );
    void                removeElement( const value_type & elem );

    /**
     *  Returns whether an element is a member of this subset
     */
    bool                operator[]( const value_type & elem );

    /**
     *  Returns whether bit index is set; equivalent to determining whether
     *  a known element is a member of this subset
     */
    bool                check_state( index_type idx ) const;

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
     *  block_iterator for subset
     */
    block_iterator begin();
    block_iterator end();

    cblock_iterator begin() const;
    cblock_iterator end() const;

    template < class E, class B, class BM, class EK >
    friend bool operator==(const variable_subset< E, B, BM, EK> & lhs, const variable_subset< E, B, BM, EK> & rhs );

    template < class E, class B, class BM, class EK >
    friend std::ostream & operator<<( std::ostream & lhs, const variable_subset< E, B, BM, EK> & rhs );

    virtual ~variable_subset();

    CLOTHO_PROTECTED

    variable_subset(powerset_type * p);
    variable_subset(powerset_type * p, const bitset_type & b);
    variable_subset( const variable_subset & vs );

    powerset_type *       m_parent;
    bitset_type         m_data;
//    unsigned int        m_ref_count;
};

#define TEMPLATE_HEADER template < class Element, class Block, class BlockMap, class ElementKeyer >
#define SUBSET_SPECIALIZATION variable_subset< Element, Block, BlockMap, ElementKeyer >
#define POWERSET_SPECIALIZATION SUBSET_SPECIALIZATION::powerset_type

TEMPLATE_HEADER
SUBSET_SPECIALIZATION::variable_subset( powerset_type * p ) :
    m_parent(p)
//    , m_ref_count(1)
{}

TEMPLATE_HEADER
SUBSET_SPECIALIZATION::variable_subset( powerset_type * p, const bitset_type & b ) :
    m_parent(p)
    , m_data( b )
//    , m_ref_count(1)
{}

TEMPLATE_HEADER
SUBSET_SPECIALIZATION::variable_subset( const variable_subset & vs ) :
    m_parent(vs.m_parent)
    , m_data( vs.m_data )
//    , m_ref_count(1)
{}

TEMPLATE_HEADER
SUBSET_SPECIALIZATION::~variable_subset() {}

TEMPLATE_HEADER
typename SUBSET_SPECIALIZATION::powerset_type * SUBSET_SPECIALIZATION::getParent() const {
    return m_parent;
}

TEMPLATE_HEADER
bool SUBSET_SPECIALIZATION::isSameFamily( pointer p ) const {
    return this->m_parent == p->m_parent;
}

TEMPLATE_HEADER
typename SUBSET_SPECIALIZATION::pointer SUBSET_SPECIALIZATION::clone() const {
    pointer sub = this->m_parent->create_subset( m_data );
    return sub;
}

//TEMPLATE_HEADER
//void SUBSET_SPECIALIZATION::copy() {
//    ++m_ref_count;
//}
//
//TEMPLATE_HEADER
//void SUBSET_SPECIALIZATION::release() {
//    --m_ref_count;
//}
//
//TEMPLATE_HEADER
//unsigned int SUBSET_SPECIALIZATION::ref_count() const {
//    return m_ref_count;
//}


TEMPLATE_HEADER
bool    SUBSET_SPECIALIZATION::operator[]( const value_type & elem ) {
    typename powerset_type::element_index_type idx = m_parent->find( elem );

    if( 0 <=  idx && idx < m_data.size() ) {
        return m_data[idx];
    }
    return false;
}

TEMPLATE_HEADER
void    SUBSET_SPECIALIZATION::addElement( const value_type & elem ) {
    std::pair< typename powerset_type::element_index_type, bool > res = m_parent->find_or_create( elem );

    if( res.first == powerset_type::npos ) return;

    if( res.first >= m_data.size() ) {
        m_data.resize( res.first + 1, false );
    }
    m_data[ res.first ] = true;
}

TEMPLATE_HEADER
void SUBSET_SPECIALIZATION::removeElement( const value_type & elem ) {
    typename powerset_type::element_index_type idx = m_parent->find( elem );

    if( idx == powerset_type::npos || m_data.size() <= idx ) return;

    m_data[idx] = false;
}

TEMPLATE_HEADER
size_t SUBSET_SPECIALIZATION::count( ) const {
    return m_data.count();
}

TEMPLATE_HEADER
size_t SUBSET_SPECIALIZATION::num_blocks( ) const {
    return m_data.num_blocks();
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
bool SUBSET_SPECIALIZATION::check_state( index_type idx ) const {
    return m_data[idx];
}

TEMPLATE_HEADER
typename SUBSET_SPECIALIZATION::block_iterator SUBSET_SPECIALIZATION::begin() {
    return m_data.m_bits.begin();
}

TEMPLATE_HEADER
typename SUBSET_SPECIALIZATION::block_iterator SUBSET_SPECIALIZATION::end() {
    return m_data.m_bits.end();
}

TEMPLATE_HEADER
typename SUBSET_SPECIALIZATION::cblock_iterator SUBSET_SPECIALIZATION::begin() const {
    return m_data.m_bits.begin();
}

TEMPLATE_HEADER
typename SUBSET_SPECIALIZATION::cblock_iterator SUBSET_SPECIALIZATION::end() const {
    return m_data.m_bits.end();
}

/**
 * Equality of variable_subsets.
 *
 * Note the dynamic_bitset::operator== compares size() of bitsets; if they are not equal the
 * sets.  However, we only attempt to maintain minimal subsets.  That is, we allow a variable
 * length zero extention of a subset to exist.  Therefore, we compare the
 * common prefix length of blocks, then require the longer sequence to be all 0s there after.
 */
TEMPLATE_HEADER
inline bool operator==( const SUBSET_SPECIALIZATION & lhs, const SUBSET_SPECIALIZATION & rhs ) {
    if( lhs.m_parent != rhs.m_parent ) return false;

    typename SUBSET_SPECIALIZATION::cblock_iterator lit = lhs.m_data.m_bits.begin(), lend = lhs.m_data.m_bits.end();
    typename SUBSET_SPECIALIZATION::cblock_iterator rit = rhs.m_data.m_bits.begin(), rend = rhs.m_data.m_bits.end();

    bool is_eq = true;
    while( is_eq ) {
        if( lit == lend ) {
            while( is_eq && rit != rend) {
                is_eq = ((*rit++) == 0);
            }
            break;
        }

        if( rit == rend ) {
            while(is_eq && lit != lend ) {
                is_eq = ((*lit++) == 0);
            }
            break;
        }

        is_eq = ((*lit++) == (*rit++));
    }
    return is_eq;
}

TEMPLATE_HEADER
std::ostream & operator<<( std::ostream & lhs, const SUBSET_SPECIALIZATION & rhs ) {
    lhs << "{" << rhs.m_ref_count << ";" << rhs.m_data << "}";
    return lhs;
}

}   // namespace powersets
}   // namespace clotho
#endif  // VARIABLE_SUBSET_HPP_
