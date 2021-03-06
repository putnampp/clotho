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
