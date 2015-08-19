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
#include <cmath>

//#include "clotho/utility/log_helper.hpp"
//#include <boost/property_tree/json_parser.hpp>
//#include <iostream>

namespace clotho {
namespace utility {

template < class URNG/*, class Result, class Tag*/ >
class random_generator< URNG, clotho::classifiers::region_classifier< qtl_allele/*, Result, Tag*/ > > {
public:
    typedef clotho::classifiers::region_classifier< qtl_allele/*, Result, Tag*/ >   result_type;
    typedef random_generator< URNG, result_type > self_type;

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
        unsigned int n = ((m_bSkip) ? 0 : m_dist(*m_rng) );

        generate_ordered_list( p, n );

//        boost::property_tree::ptree l;
//        add_value_array(l, p.begin(), p.end() );
//        boost::property_tree::write_json( std::cerr, l );

        return result_type(p);
    }

protected:

    void method1(typename result_type::param_type & p, unsigned int n) {
        if( n == 0 ) return;

        std::set< double > positions;
        while( positions.size() < n ) {
            positions.insert( m_uniform( *m_rng ) );
        }

        qtl_allele::trait_weights coeff;
        for( std::set< double >::iterator it = positions.begin(); it != positions.end(); ++it ) {
            p.push_back( qtl_allele( *it, DEFAULT_SELECTION, DEFAULT_DOMINANCE, DEFAULT_NEUTRAL, 0, coeff) );
        }
    }

    // Linear generation of ordered random allele locations
    // @article{Bentley:1980:GSL:355900.355907,
    //  author = {Bentley, Jon Louis and Saxe, James B.},
    //  title = {Generating Sorted Lists of Random Numbers},
    //  journal = {ACM Trans. Math. Softw.},
    //  issue_date = {Sept. 1980},
    //  volume = {6},
    //  number = {3},
    //  month = sep,
    //  year = {1980},
    //  issn = {0098-3500},
    //  pages = {359--364},
    //  numpages = {6},
    //  url = {http://doi.acm.org/10.1145/355900.355907},
    //  doi = {10.1145/355900.355907},
    //  acmid = {355907},
    //  publisher = {ACM},
    //  address = {New York, NY, USA},
    //  }
    void generate_ordered_list( typename result_type::param_type & p, unsigned int N ) {
        double curmax = 0.0;
        qtl_allele::trait_weights coeff;
        while( N ) {
            // generates list in reverse order
            // hence taking the complement
            curmax += (log( m_uniform( *m_rng ) ) / (double)(N));
            p.push_back( qtl_allele( (1.0 - exp(curmax)), DEFAULT_SELECTION, DEFAULT_DOMINANCE, DEFAULT_NEUTRAL, 0, coeff) );
            --N;
        }
    }

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
