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
#ifndef BINARY_CLASSIFIER_HPP_
#define BINARY_CLASSIFIER_HPP_

namespace clotho {
namespace classifiers {

template < class Classifier, class Tag >
class binary_classifier {
public:
    typedef Classifier  classifier_type;
    typedef typename classifier_type::element_type  element_type;
    typedef bool        result_type;

    typedef binary_classifier< Classifier, Tag > self_type;

    binary_classifier( ) {}

    binary_classifier( const classifier_type & c ) : m_cfier( c ) {}

    binary_classifier( const self_type & o ) : m_cfier( o.m_cfier ) {}

    result_type operator()( element_type const& elem ) {
        typename classifier_type::result_type res = m_cfier( elem );

        return Tag::eval( res );
    }

    template < class Iter >
    result_type operator()( Iter first, size_t offset ) {
        return operator()( *(first + offset) );
    }

    size_t event_count() const {
        return m_cfier.event_count();
    }

    virtual ~binary_classifier() {}

protected:
    classifier_type m_cfier;
};

}   // namespace classifiers
}   // namespace clotho

#include "clotho/utility/random_generator.hpp"

namespace clotho {
namespace utility {

template < class URNG, class Classifier, class Tag >
class random_generator< URNG, clotho::classifiers::binary_classifier< Classifier, Tag > > : public random_generator< URNG, Classifier > {
public:
    typedef clotho::classifiers::binary_classifier< Classifier, Tag > result_type;
    typedef random_generator< URNG, result_type >   self_type;

    typedef random_generator< URNG, Classifier >    base_type;

    random_generator( URNG & rng, boost::property_tree::ptree & config ) : base_type( rng, config ) {}

    result_type operator()() {
        return result_type( ((base_type *)this)->operator()());
    }

    virtual ~random_generator() {}
};

}   // namespace utility {
}   // namespace clotho {

#endif  // BINARY_CLASSIFIER_HPP_
