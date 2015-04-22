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
#ifndef BERNOULLI_CLASSIFIER_HPP_
#define BERNOULLI_CLASSIFIER_HPP_

#include <boost/random/bernoulli_distribution.hpp>

namespace clotho {
namespace classifiers {


template < class RNG, class RealType = double, class Result = bool >
class bernoulli_classifier {
public:
    typedef RealType real_type;
    typedef Result  result_type;

    typedef bernoulli_classifier< RNG, RealType, Result > self_type;

    typedef boost::random::bernoulli_distribution< RealType > bernoulli_type;

    bernoulli_classifier( RNG * rng, real_type p = 0.5 ) :
        m_rng( rng )
        , m_bern( p ) {
    }

    bernoulli_classifier( const self_type & rhs ) : m_rng( rhs.m_rng ), m_bern( rhs.m_bern.param()) {}

    template < class Element >
    result_type operator()( const Element & elem ) const {
        return m_bern( *m_rng );
    }

    template < class ElementIterator >
    result_type operator()( ElementIterator elem_it, size_t idx) const {
        return m_bern( *m_rng );
    }

protected:
    RNG             * m_rng;
    bernoulli_type  m_bern;
};

}   // namespace classifier {
}   // namespace clotho {

#endif  // BERNOULLI_CLASSIFIER_HPP_
