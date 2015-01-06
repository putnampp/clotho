#ifndef CLOTHO_ELEMENT_ITERATOR_HPP_
#define CLOTHO_ELEMENT_ITERATOR_HPP_

#include <utility>

namespace clotho {
namespace utility {

template < class Container >
class element_iterator {
public:
    typedef element_iterator< Container >   self_type;
    typedef Container                       container_type;
    typedef typename Container::value_type  value_type;

    typedef typename Container::iterator    iterator;

    element_iterator( container_type & con ) :
        m_container(&con)
        , m_it( con.begin() ) {
    }

    element_iterator( const self_type & other) :
        m_container( other.m_container )
        , m_it( other.m_it ) {
    }

    element_iterator< Container > & operator++() {
        ++m_it;
    }

    element_iterator< Container > operator++( int ) {
        element_iterator< Container > tmp( *this );
        operator++();
        return tmp;
    }

    value_type & operator*() {
        return (*m_it);
    }

    bool operator==( const self_type & rhs ) const {
        return (m_it == rhs.m_it);
    }

    bool operator!=( const self_type & rhs ) const {
        return (m_it != rhs.m_it);
    }

    virtual ~element_iterator() {}

protected:
    container_type  * m_container;
    iterator        m_it;
};

enum pair_iterator_state { FIRST, SECOND, DONE };

template < class Element >
class element_iterator< std::pair< Element, Element > > {
public:
    typedef Element element_type;
    typedef std::pair< element_type, element_type > pair_type;

    typedef element_iterator< pair_type > self_type;

    element_iterator( pair_type & p, pair_iterator_state s = pair_iterator_state::DONE ) :
        m_elements( &p )
        , m_state( s ) {
    }

    element_iterator( pair_type * p, pair_iterator_state s = pair_iterator_state::DONE ) :
        m_elements( p )
        , m_state( s ) {
    }

    element_iterator( const self_type & other ) :
        m_elements( other.m_elements )
        , m_state( other.m_state ) {
    }

    self_type & operator++() {
        switch( m_state ) {
        case FIRST:
            m_state = SECOND;
            break;
        case SECOND:
            m_state = DONE;
            break;
        default:
            break;
        }
        return *this;
    }

    self_type operator++(int) {
        self_type res(*this);
        operator++();
        return res;
    }

    element_type & operator*() {
        return ((m_state == FIRST) ? m_elements->first : m_elements->second);
    }

    bool operator==( const self_type & rhs ) const {
        return m_elements == rhs.m_elements && m_state == rhs.m_state;
    }

    bool operator!=( const self_type & rhs ) const {
        return ( m_elements != rhs.m_elements || m_state != rhs.m_state);
    }

protected:
    pair_type   * m_elements;
    pair_iterator_state m_state;
};

}   // namespace utility {
}   // namespace clotho {
#endif  // ELEMENT_ITERATOR_HPP_
