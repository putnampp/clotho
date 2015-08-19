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
#ifndef QTL_HASH_CLASSIFIER_HPP_
#define QTL_HASH_CLASSIFIER_HPP_

#include "qtl_allele.h"

#include "clotho/utility/random_generator.hpp"
#include "clotho/classifiers/hash_classifier.hpp"
#include "clotho/classifiers/helpers/normalize.hpp"

#include <boost/random/poisson_distribution.hpp>

namespace clotho {
namespace classifiers {
namespace helpers {

// used in std::hash< clotho::classifiers::hash_classifier< Element, Bins>> specialization
inline double normalize_arg( qtl_allele const& q ) {
    return q.m_key;
}

}   // namespace helpers {
}   // namespace classifiers {
}   // namespace clotho {


namespace clotho {
namespace utility {

template < class URNG, unsigned int Bins >
class random_generator< URNG, clotho::classifiers::hash_classifier< qtl_allele, Bins > > {
public:
    typedef clotho::classifiers::hash_classifier< qtl_allele, Bins > result_type;

    typedef random_generator< URNG, result_type > self_type;

    typedef boost::random::uniform_01< double > position_dist_type;
    typedef boost::random::poisson_distribution< unsigned int, double > event_dist_type;

    random_generator( URNG & rng, boost::property_tree::ptree & config ) :
        m_rng( &rng )
        , m_lambda( DEFAULT_EVENT_RATE )
    {
        parseConfig( config );
    }

    result_type operator()() {
        result_type res;

        unsigned int N = m_event_dist( *m_rng );

        typename result_type::param_type::bin_type _sum = 0;
        unsigned int bin_idx = 0;
        double accum = 0.;
        qtl_allele::trait_weights coeff;

        while( N ) {
            accum += (log( m_pos_dist( *m_rng ) ) / (double) N);

            res.param()->class_bounds.push_back( qtl_allele( (1.0 - exp(accum)), DEFAULT_SELECTION, DEFAULT_DOMINANCE, DEFAULT_NEUTRAL, 0, coeff) );

            unsigned int bidx = m_hash( res.params()->class_bounds.back() );

            while( bin_idx < bidx ) {
                res.param()->bin_psum[ ++bin_idx ] = _sum;
            }

            ++_sum;
            --N;
        }

        while( bin_idx <= result_type::BINS ) {
            res.param()->bin_psum[ bin_idx++ ] = _sum;
        }

        return res;
    }

    virtual ~random_generator() {}
protected:
    URNG * m_rng;
    position_dist_type              m_pos_dist;
    event_dist_type                 m_event_dist;
    typename result_type::hash_type m_hash;
};

}   // namespace utility
}   // namespace clotho

#endif  // QTL_HASH_CLASSIFIER_HPP_
