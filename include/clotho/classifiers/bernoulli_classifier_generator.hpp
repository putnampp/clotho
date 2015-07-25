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
#ifndef BERNOULLI_CLASSIFIER_GENERATOR_HPP_
#define BERNOULLI_CLASSIFIER_GENERATOR_HPP_

#include "clotho/utility/random_generator.hpp"
#include "clotho/classifiers/bernoulli_classifier.hpp"

namespace clotho {
namespace utility {

template < class URNG >
class random_generator< URNG, clotho::classifiers::bernoulli_classifier< URNG, double, bool > > {
public:
    typedef random_generator< URNG, clotho::classifiers::bernoulli_classifier< URNG, double, bool > > self_type;
    typedef clotho::classifiers::bernoulli_classifier< URNG, double, bool >   result_type;

    random_generator( URNG & rng, boost::property_tree::ptree & config ) :
        m_rng( &rng )
        , m_p( 0.5 ) {
        parseConfig( config );
    }

    random_generator( URNG & rng, double p = 0.5 ) :
        m_rng( &rng )
        , m_p( p ) {
    }

    result_type operator()() {
        return result_type( m_rng, m_p );
    }

protected:

    void parseConfig( boost::property_tree::ptree & config ) {
        const string ppath = "classifiers.toolkit.bernoulli.p";

        if( config.get_child_optional( ppath ) == boost::none ) {
            config.put( ppath, m_p );
        } else {
            m_p = config.get< double >( ppath );

            assert( 0. <= m_p && m_p <= 1.);
        }
    }

    URNG        * m_rng;
    double      m_p;
};

}   // namespace utility {
}   // namespace clotho {

#endif  // BERNOULLI_CLASSIFIER_GENERATOR_HPP_
