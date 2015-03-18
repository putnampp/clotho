#ifndef QTL_ALLELE_GENERATOR_HPP_
#define QTL_ALLELE_GENERATOR_HPP_

#include "qtl_allele.h"
#include "clotho/utility/random_generator.hpp"

#include <boost/random/uniform_real_distribution.hpp>
#include <algorithm>

namespace clotho {
namespace utility {

template < class URNG >
class random_generator< URNG, qtl_allele > {
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
        m_rng( &rng )
        , m_alleles(rng, config )
        , m_traits( rng, config )
        , m_nTraits( 1 ) {
        parseConfig( config );
    }

    random_generator( URNG & rng, unsigned int n = 1, double trait_mean = 0.0, double trait_sigma = 1.0 ) :
        m_rng( &rng )
        , m_alleles( rng )
        , m_traits(rng, trait_mean, trait_sigma)
        , m_nTraits(n) {
    }

    result_type operator()( unsigned int age = 0 ) {
        basic_allele all = m_alleles( age );
        typename qtl_allele::trait_weights W;

        std::generate_n( std::back_inserter(W), m_nTraits, m_traits);

        return qtl_allele(all, W );
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

    URNG            * m_rng;
    random_allele   m_alleles;
    random_traits   m_traits;

    unsigned int    m_nTraits;
};

}   // namespace utility {
}   // namespace clotho {

#endif  // QTL_ALLELE_GENERATOR_HPP_
