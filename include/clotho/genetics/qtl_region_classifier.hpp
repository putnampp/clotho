#ifndef QTL_REGION_CLASSIFIER_HPP_
#define QTL_REGION_CLASSIFIER_HPP_

#include "qtl_allele.h"
#include "clotho/utility/random_generator.hpp"
#include "clotho/classifiers/region_classifier.hpp"

#include <boost/random/poisson_distribution.hpp>

namespace clotho {
namespace utility {

template < class URNG, class Result, class Tag >
class random_generator< URNG, clotho::classifiers::region_classifier< qtl_allele, Result, Tag > > {
public:
    typedef random_generator< URNG, clotho::classifiers::region_classifier< qtl_allele, Result, Tag > > self_type;
    typedef clotho::classifiers::region_classifier< qtl_allele, Result, Tag >   result_type;

    typedef boost::random::uniform_01< double >                                 uniform_type;
    typedef boost::random::poisson_distribution< unsigned int, double >         dist_type;

    random_generator( URNG & rng, boost::property_tree::ptree & config ) :
        m_rng( &rng )
        , m_dist( DEFAULT_RECOMB_RATE ) {
        parseConfig( config );
    }

    random_generator( URNG & rng, double mu = DEFAULT_RECOMB_RATE ) :
        m_rng( &rng ), m_dist( mu ) {
    }

    result_type operator()() {
        typename result_type::param_type p;

        unsigned int n = m_dist( *m_rng );

        for( unsigned int i = 0; i < n; ++i ) {
            double k = m_uniform( *m_rng );
            qtl_allele::trait_weights coeff;

            p.push_back( qtl_allele( k, DEFAULT_SELECTION, DEFAULT_DOMINANCE, DEFAULT_NEUTRAL, coeff) );
        }

        std::sort( p.begin(), p.end() );
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
            typename dist_type::param_type p( m );

            m_dist.param( p );
        }
    }

    URNG        * m_rng;
    uniform_type    m_uniform;
    dist_type       m_dist;
};

}   // namespace utility {
}   // namespace clotho {

#endif  // QTL_REGION_CLASSIFIER_HPP_
