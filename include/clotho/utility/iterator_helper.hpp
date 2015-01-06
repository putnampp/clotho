#ifndef CLOTHO_ITERATOR_HELPER_HPP_
#define CLOTHO_ITERATOR_HELPER_HPP_

#include "clotho/utility/element_iterator.hpp"

namespace clotho {
namespace utility {

template < class Container >
struct iterator_helper {
    typedef clotho::utility::element_iterator< Container >  type;
    typedef std::pair< type, type >                         range_type;

    static type make_first( Container & con ) {
        return type( con );
    }

    static type make_last( Container & con ) {
        return type( con );
    }

    static range_type make_range( Container & con ) {
        return std::make_pair( type(con), type(con) );
    }
};

template < class Element >
struct iterator_helper< std::pair< Element, Element > > {
    typedef clotho::utility::element_iterator< std::pair< Element, Element > > type;
    typedef std::pair< type, type >                                            range_type;

    static type    make_first( std::pair< Element, Element > & ind ) {
        return type( ind, pair_iterator_state::FIRST);
    }

    static type    make_last( std::pair< Element, Element > & ind ) {
        return type( ind );
    }

    static range_type make_range( std::pair< Element, Element > & ind ) {
        return std::make_pair( type(ind, pair_iterator_state::FIRST), type(ind) );
    }
};

template < class Element >
struct iterator_helper< std::shared_ptr< Element > > :
        public iterator_helper< Element > {
    typedef iterator_helper< Element >      base_type;
    typedef typename base_type::type        type;
    typedef typename base_type::range_type  range_type;

    static type make_first( std::shared_ptr< Element > s ) {
        return base_type::make_first( *s );
    }

    static type make_last( std::shared_ptr< Element > s ) {
        return base_type::make_last( *s );
    }

    static range_type make_range( std::shared_ptr< Element > s ) {
        return base_type::make_range( *s );
    }
};

}   // namespace utility {
}   // namespace clotho {

#endif  // CLOTHO_ITERATOR_HELPER_HPP_
