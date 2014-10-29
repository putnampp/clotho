#ifndef CLOTHO_ELEMENT_GENERATOR_HPP_
#define CLOTHO_ELEMENT_GENERATOR_HPP_

#include <cassert>

namespace clotho {
namespace mutations {

template < class Element, class Generator >
class element_generator {
public:
    typedef Element                         element_type;
    typedef Generator                       generator_type;

    element_generator( generator_type & g ) :
        m_gen( g )
    {}

    element_type operator()( ) {
        element_type elem(m_gen() );
        return elem;
    }

    template < class URNG >
    element_type operator()( URNG & r ) {
        element_type elem( m_gen(r) );
        return elem;
    }

    template < class URNG, class Pred >
    element_type operator()( URNG & r, Pred * p ) {
        element_type elem( m_gen(r) );

        while( !(*p)( elem ) ) {
            elem = element_type( m_gen(r) );
        }

        return elem;
    }

    template < class URNG, class ResultIterator >
    void operator()( URNG & r, unsigned int nMut, ResultIterator res_it ) {
        while( nMut-- ) {
            (*res_it++) = element_type( m_gen(r) );
        }
    }

    template < class URNG, class Pred, class ResultIterator >
    void operator()( URNG & r, unsigned int nMut, Pred * p, ResultIterator res_it ) {
        while( nMut-- ) {
            element_type elem( m_gen(r) );
            if( !(*p)( elem ) ) {
                ++nMut;
            } else {
                (*res_it++) = elem;
            }
        }
    }

    virtual ~element_generator() {}
protected:
    generator_type  m_gen;
};

}   // namespace mutations {
}   // namespace clotho {

#endif  // CLOTHO_ELEMENT_GENERATOR_HPP_