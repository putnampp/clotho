#ifndef VARIABLE_SUBSET_HPP_
#define VARIABLE_SUBSET_HPP_

#include "clotho/powersets/powerset.hpp"
#include "boost/dynamic_bitset.hpp"

namespace clotho {
namespace powersets {

template < class Element, class Block = unsigned long, class BlockMap = block_map< Element, Block>, class ElementKeyer = element_key_of< Element > >
class variable_subset {
    public:
        typedef Element value_type;

        typedef clotho::powersets::powerset< Element, Block, BlockMap, ElementKeyer, variable_subset< Element, Block, BlockMap, ElementKeyer > > powerset_type;
        typedef boost::dynamic_bitset< unsigned long > bitset_type;

        friend class clotho::powersets::powerset< Element, Block, BlockMap, ElementKeyer, variable_subset< Element, Block, BlockMap, ElementKeyers > >;

        variable_subset *   clone() const;
        variable_subset *   copy();

        void                release();
        unsigned int        ref_count() const;

        void                addElement( const Element & elem );
        void                removeElement( const Element & elem );

        friend bool operator==(const variable_subset & lhs, const variable_subset & rhs );

        virtual ~variable_subset();
    protected:
        variable_subset(powerset_type * p);
        variable_subset( const variable_subset & vs );

        powerset_type *       m_parent;
        bitset_type         m_data;
        unsigned int        m_ref_count;
    };

#define TEMPLATE_HEADER template < class Element, class Block, class BlockMap, class ElementKeyer >
#define SUBSET_SPECIALIZATION variable_subset< Element, Block, BlockMap, ElementKeyer >
#define POWERSET_SPECIALIZATION SUBSET_SPECIALIZATION::powerset_type

TEMPLATE_HEADER
SUBSET_SPECIALIZATION::variable_subset( powerset_type * p ) :
    m_parent(p)
    , m_ref_count(1)
{}

TEMPLATE_HEADER
SUBSET_SPECIALIZATION::variable_subset( const variable_subset & vs ) :
    m_parent(vs.m_parent)
    , m_data( vs.m_data )
    , m_ref_count(1)
{}

TEMPLATE_HEADER
SUBSET_SPECIALIZATION::~variable_subset() {}

TEMPLATE_HEADER
typename SUBSET_SPECIALIZATION * SUBSET_SPECIALIZATION::clone() const {
    variable_subset * sub = this->m_parent->clone_subset( this );

    return sub;
}

TEMPLATE_HEADER
typename SUBSET_SPECIALIZATION * SUBSET_SPECIALIZATION::copy() {
    m_parent->copy_subset( this );
    return this;
}

TEMPLATE_HEADER
void SUBSET_SPECIALIZATION::release() {
    m_parent->release_subset( this );
}

TEMPLATE_HEADER
void    SUBSET_SPECIALIZATION::addElement( const POWERSET_SPECIALIZATION::value_type & elem ) {
    typename powerset_type::element_index_type idx = m_parent->find_or_create( elem );

    if( idx == powerset_type::npos ) return;

    if( idx >= m_data.size() ) {
        m_data.resize( idx + 1, false );
    }
    m_data[ idx ] = true;
}

TEMPLATE_HEADER
void SUBSET_SPECIALIZATION::removeElement( const value_type & elem ) {
    typename powerset_type::element_index_type idx = m_parent->find( elem );

    if( idx == powerset_type::npos || m_data.size() <= idx ) return;

    m_data[idx] = false;
}

}   // namespace powersets
}   // namespace clotho
#endif  // VARIABLE_SUBSET_HPP_
