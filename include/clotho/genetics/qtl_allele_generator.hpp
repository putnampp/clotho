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
#ifndef QTL_ALLELE_GENERATOR_HPP_
#define QTL_ALLELE_GENERATOR_HPP_

#include "qtl_allele.h"
#include "clotho/utility/random_generator.hpp"

#include <boost/random/uniform_real_distribution.hpp>
#include <algorithm>

namespace clotho {
namespace utility {

template < class URNG >
class random_generator< URNG, qtl_allele > : public clotho::utility::random_generator< URNG, basic_allele > {
public:
    typedef random_generator< URNG, qtl_allele >    self_type;
    typedef qtl_allele                              result_type;

    typedef double real_type;

    typedef real_type key_type;
    typedef real_type selection_type;
    typedef real_type dominance_type;

    typedef boost::random::uniform_01< real_type > dist_type;

    typedef clotho::utility::random_generator< URNG, basic_allele >                  random_allele;
    typedef clotho::utility::random_generator< URNG, qtl_allele::weight_type >       random_traits;
    typedef typename random_traits::param_type  param_type;

    random_generator( URNG & rng, boost::property_tree::ptree & config ) :
//        m_rng( &rng )
//        , m_alleles(rng, config )
        random_allele(rng, config )
        , m_traits( rng, config )
        , m_nTraits( 1 ) {
        parseConfig( config );
    }

    random_generator( URNG & rng, unsigned int n = 1, double trait_mean = 0.0, double trait_sigma = 1.0 ) :
//        m_rng( &rng )
//        , m_alleles( rng )
        random_allele( rng )
        , m_traits(rng, trait_mean, trait_sigma)
        , m_nTraits(n) {
    }

    result_type operator()( unsigned int age = 0 ) {
//        basic_allele all = m_alleles( age );
//        typename qtl_allele::trait_weights W;
//
//        qtl_allele q(all, W );
//
        qtl_allele q;

        random_allele::generate(q, age);

        if( q.isNeutral() ) {
            // neutral alleles, by definition, do not contribute
            // to traits (phenotype)
            std::fill_n( std::back_inserter(q.m_weights), m_nTraits, 0);
        } else {
            std::generate_n( std::back_inserter(q.m_weights), m_nTraits, m_traits);
        }

        return q;
    }

protected:

    void parseConfig( boost::property_tree::ptree & config ) {
        std::ostringstream oss;
        oss /*<< CONFIG_BLOCK_K << "."*/ << TRAIT_BLOCK_K << "." << MAX_K;

        if( config.get_child_optional( oss.str() ) == boost::none ) {
            config.put( oss.str(), m_nTraits );
        } else {
            m_nTraits = config.get< unsigned int >( oss.str(), 1 );
        }
    }

//    URNG            * m_rng;
//    random_allele   m_alleles;
    random_traits   m_traits;

    unsigned int    m_nTraits;
};

}   // namespace utility {
}   // namespace clotho {

#endif  // QTL_ALLELE_GENERATOR_HPP_
