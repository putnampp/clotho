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

//#include "clotho/utility/log_helper.hpp"
//#include <boost/property_tree/json_parser.hpp>
//#include <iostream>

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
        unsigned int n = ((m_bSkip) ? 0 : m_dist(*m_rng) );

        generate_ordered_list( p, 0.0, 1.0, n );

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

/**
 * Construct a position ordered list of qtl_alleles in linear time
 * Avoids potential of duplicate positions in the list
 * Avoids having to sort the list
 */
//   [lo, hi)        | n   | r    | val                        | lside           | rside | rec | p
//   [0, 1)          | 10  | 0.5  | 0.5 * 1 + 0 = 0.5          | (0.5 * 9) = 4   | 5     | l   |
//   [0, 0.5)        | 4   | 0.3  | 0.3 * (0.5) + 0 = 0.15     | (0.3 * 3) = 0   | 3     | l   |
//   [0, 0.15)       | 0   | -    | -                          | -               | -     | -   |
//   [0, 0.5)        | 3   | 0.3  | 0.3 * (0.5) + 0 = 0.15     | (0.3 * 2) = 0   | 2     | r   | 0.15
//   [0.15, 0.5 )    | 2   | 0.8  | 0.8 * (0.35) +0.15 = 0.43  | (0.8 * 1) = 0   | 1     | l   | 0.15
//   [0.15, 0.43)    | 0   | -    | -                          | -               | -     | -   | 0.15
//   [0.43, 0.5)     | 1   | 0.11 | 0.11* (0.43) +0.15 = 0.4377| 0               | 0     | r   | 0.15, 0.43
//   [0.43, 0.4377)  | 0   | -    | -                          | -               | -     | -   | 0.15, 0.43
//   [0.43, 0.5)     | 1   | 0.11 | 0.11* (0.43) +0.15 = 0.4377| 0               | 0     | r   | 0.15, 0.43, 0.4377
//   [0.4377, 0.5)   | 0   | -    | -                          | -               | -     | -   | 0.15, 0.43, 0.4377
//   [0, 1)          | 10  | 0.5  | 0.5 * 1 + 0 = 0.5          | (0.5 * 9) = 4   | 5     | r   | 0.15, 0.43, 0.4377, 0.5
//   [0.5, 1)        | 5   | 0.9  | 0.9 * (0.5) + 0.5 = 0.95   | (0.5 * 4) = 2   | 2     | l   | 0.15, 0.43, 0.4377, 0.5
//   [0.5, 0.95)     | 2   | 0.02 | 0.02 * (0.45) + 0.5 = 0.509| (0.02 * 1 ) = 0 | 1     | l   | 0.15, 0.43, 0.4377, 0.5
//   [0.5, 0.509)    | 0   | -    | -                          | -               | -     | -   | 0.15, 0.43, 0.4377, 0.5
//   [0.5, 0.95)     | 2   | 0.02 | 0.02 * (0.45) + 0.5 = 0.509| (0.02 * 1 ) = 0 | 1     | r   | 0.15, 0.43, 0.4377, 0.5, 0.509

    void generate_ordered_list( typename result_type::param_type & p, double lo, double hi, unsigned int n ) {
        if( n-- == 0 ) return;
        
        // m_uniform produces X ~ [0, 1)
        double r = m_uniform( *m_rng );

        while( r == 0.0 || r == 1.0 ) { r = m_uniform( *m_rng ); }

        double val = r * (hi - lo) + lo;

        unsigned int lside = r * n;

        generate_ordered_list( p, lo, val, lside);

        qtl_allele::trait_weights coeff;
        p.push_back( qtl_allele( val, DEFAULT_SELECTION, DEFAULT_DOMINANCE, DEFAULT_NEUTRAL, 0, coeff) );
        
        generate_ordered_list( p, val, hi, n - lside );
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
