#ifndef CLOTHO_ITERATOR_HELPER_HPP_
#define CLOTHO_ITERATOR_HELPER_HPP_

#include "clotho/utility/element_iterator.hpp"

namespace clotho {
namespace utility {

template < class Container >
struct iterator_helper {
    typedef clotho::utility::element_iterator< Container > type;

    static type make_first( Container & con ) {
        return type( con );
    }

    static type make_second( Container & con ) {
        return type( con );
    }
};

template < class Element >
struct iterator_helper< std::pair< Element, Element > > {
    typedef element_iterator< std::pair< Element, Element > > type;

    static type    make_first( std::pair< Element, Element > & ind ) {
        return element_iterator< Element >( ind, pair_iterator_state::FIRST);
    }

    static type    make_last( std::pair< Element, Element > & ind ) {
        return element_iterator< Element >( ind );
    }
};

#include <memory>

template < class Element >
struct iterator_helper< std::shared_ptr< Element > > {
    typedef typename iterator_helper< Element >::type type;

    static type make_first( std::shared_ptr< Element > s ) {
        return iterator_helper< Element >::make_first( *s );
    }

    static type make_last( std::shared_ptr< Element > s ) {
        return iterator_helper< Element >::make_last( *s );
    }
};

}   // namespace utility {
}   // namespace clotho {

#endif  // CLOTHO_ITERATOR_HELPER_HPP_
