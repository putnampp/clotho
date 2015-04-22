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
#ifndef QTL_REGION_CLASSIFIER_HPP_
#define QTL_REGION_CLASSIFIER_HPP_

#include "qtl_allele.h"
#include "clotho/utility/random_generator.hpp"
#include "clotho/classifiers/region_classifier.hpp"

#include <boost/random/poisson_distribution.hpp>
#include <set>

namespace clotho {
namespace utility {

template < class URNG, class Result, class Tag >
class random_generator< URNG, clotho::classifiers::region_classifier< qtl_allele, Result, Tag > > {
public:
    typedef random_generator< URNG, clotho::classifiers::region_classifier< qtl_allele, Result, Tag > > self_type;
    typedef clotho::classifiers::region_classifier< qtl_allele, Result, Tag >   result_type;

    typedef boost::random::uniform_01< double >                                 uniform_type;   // key distribution
    typedef boost::random::poisson_distribution< unsigned int, double >         dist_type;      // region distribution

    random_generator( URNG & rng, boost::property_tree::ptree & config ) :
        m_rng( &rng )
        , m_dist( DEFAULT_RECOMB_RATE )
        , m_bSkip(false) {
        parseConfig( config );
    }

    random_generator( URNG & rng, double mu = DEFAULT_RECOMB_RATE ) :
        m_rng( &rng ), m_dist( mu ) {
    }

    result_type operator()() {
        typename result_type::param_type p;

        unsigned int n = ((m_bSkip)? 0 : m_dist( *m_rng ));

        if( n ) {
            std::set< double > positions;
            while( positions.size() < n ) {
                positions.insert( m_uniform( *m_rng ) );
            }

            qtl_allele::trait_weights coeff;
            for( std::set< double >::iterator it = positions.begin(); it != positions.end(); ++it ) {
                p.push_back( qtl_allele( *it, DEFAULT_SELECTION, DEFAULT_DOMINANCE, DEFAULT_NEUTRAL, 0, coeff) );
            }
        }
        return result_type( p );
    }

protected:

    void parseConfig( boost::property_tree::ptree & config ) {
        std::ostringstream oss;
        oss /*<< CONFIG_BLOCK_K << "."*/ << REC_BLOCK_K << "." << RATE_PER_REGION_K;

        if( config.get_child_optional( oss.str() ) == boost::none ) {
            config.put( oss.str(), m_dist.mean() );
        } else {
            double m = config.get< double >( oss.str() );

            m_bSkip = (m == 0.0);

            if( m_bSkip ) {
                m = 0.00000000001;
            } else if( m < 0.0 ) {
                m = std::abs( m );
            }

            typename dist_type::param_type p( m );

            m_dist.param( p );
        }
    }

    URNG        * m_rng;
    uniform_type    m_uniform;
    dist_type       m_dist;
    bool        m_bSkip;
};

}   // namespace utility {
}   // namespace clotho {

#endif  // QTL_REGION_CLASSIFIER_HPP_
