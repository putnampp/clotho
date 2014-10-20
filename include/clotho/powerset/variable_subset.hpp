#ifndef VARIABLE_SUBSET_HPP_
#define VARIABLE_SUBSET_HPP_

//#include <boost/dynamic_bitset.hpp>

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

        friend class clotho::powersets::powerset< Element, variable_subset< Element, Block, BlockMap, ElementKeyer >, Block, BlockMap, ElementKeyer >;

        //self_type *   clone() const;
//        self_type *   copy();

//        void                release();
//        unsigned int        ref_count() const;

        void                addElement( const value_type & elem );
        void                removeElement( const value_type & elem );

        bool                check_state( index_type idx ) const;

        size_t              count() const;

        template < class E, class B, class BM, class EK >
        friend bool operator==(const variable_subset< E, B, BM, EK> & lhs, const variable_subset< E, B, BM, EK> & rhs );

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

//TEMPLATE_HEADER
//SUBSET_SPECIALIZATION * SUBSET_SPECIALIZATION::clone() const {
//    variable_subset * sub = this->m_parent->clone_subset( this );
//
//    return sub;
//}

//TEMPLATE_HEADER
//SUBSET_SPECIALIZATION * SUBSET_SPECIALIZATION::copy() {
//    m_parent->copy_subset( this );
//    return this;
//}

//TEMPLATE_HEADER
//void SUBSET_SPECIALIZATION::release() {
//    m_parent->release_subset( this );
//}

TEMPLATE_HEADER
void    SUBSET_SPECIALIZATION::addElement( const value_type & elem ) {
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

TEMPLATE_HEADER
size_t SUBSET_SPECIALIZATION::count( ) const {
    return m_data.count();
}

TEMPLATE_HEADER
bool SUBSET_SPECIALIZATION::check_state( index_type idx ) const {
    return m_data[idx];
}

TEMPLATE_HEADER
inline bool operator==( const SUBSET_SPECIALIZATION & lhs, const SUBSET_SPECIALIZATION & rhs ) {
    return (lhs.m_parent == rhs.m_parent && lhs.m_data == rhs.m_data);
}

}   // namespace powersets
}   // namespace clotho
#endif  // VARIABLE_SUBSET_HPP_
