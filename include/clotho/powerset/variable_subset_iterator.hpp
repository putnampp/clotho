#ifndef CLOTHO_VARIABLE_SUBSET_ITERATOR_HELPER_HPP_
#define CLOTHO_VARIABLE_SUBSET_ITERATOR_HELPER_HPP_

#include "clotho/utility/iterator_helper.hpp"
#include "clotho/powersets/variable_subset.hpp"

namespace clotho {
namespace utility {

template < class E, class B, class BK, class EK >
class element_iterator< clotho::powersets::variable_subset< E, B, BK, EK > > {
public:
    typedef clotho::powersets::variable_subset< E, B, BK, EK > container_type;
    typedef element_iterator< container_type > self_type;
    typedef E       value_type;

    typedef typename container_type::bitset_type::size_type size_type;

    typedef typename container_type::powerset_type::variable_iterator  iterator;

    element_iterator( container_type & con, size_type idx = container_type::bitset_type::npos ) :
        m_container( &con )
        , m_idx(idx)
        , m_var(con.getParent()->variable_begin())
    { }

    element_iterator( const self_type & other ) :
        m_container( other.m_container )
        , m_idx( other.m_idx )
        , m_var( other.m_var )
    {}

    self_type & operator++() {
        if( m_idx != container_type::bitset_type::npos ) {
            m_idx = m_container->find_next(m_idx);
        }
        return *this;
    }

    self_type operator++( int ) {
        self_type tmp(*this);
        operator++();
        return tmp;
    }

    value_type &  operator*() {
        assert( m_idx != container_type::bitset_type::npos );

        return (*(m_var + m_idx));
    }

    bool operator==( const self_type & rhs ) {
        return (m_container == rhs.m_container && m_idx == rhs.m_idx);
    }

    bool operator!=( const self_type & rhs ) {
        return (m_container != rhs.m_container || m_idx != rhs.m_idx);
    }

    virtual ~element_iterator() {}
protected:
    container_type  * m_container;
    size_type       m_idx;
    iterator        m_var;
};

template < class E, class B, class BK, class EK >
struct iterator_helper< clotho::powersets::variable_subset< E, B, BK, EK > > {
    typedef clotho::powersets::variable_subset< E, B, BK, EK > subset_type;
    typedef element_iterator< subset_type > type;

    static type make_first( subset_type & s ) {
        return type( s, s.find_first() );
    }

    static type make_last( subset_type & s ) {
        return type(s);
    }
};

}   // namespace utility {
}   // namespace clotho {
#endif  // CLOTHO_VARIABLE_SUBSET_ITERATOR_HELPER_HPP_
