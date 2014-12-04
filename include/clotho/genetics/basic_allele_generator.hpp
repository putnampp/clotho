#ifndef BASIC_ALLELE_GENERATOR_HPP_
#define BASIC_ALLELE_GENERATOR_HPP_

#include "basic_allele.h"
#include "clotho/utility/random_generator.hpp"

#include <boost/random/uniform_01.hpp>
#include <boost/random/bernoulli_distribution.hpp>

extern const string NEUTRAL_P_K;

namespace clotho {
namespace utility {

template < class URNG >
class random_generator< URNG, basic_allele > {
public:
    typedef random_generator< URNG, basic_allele >  self_type;
    typedef basic_allele                            result_type;

    typedef double real_type;

    typedef real_type key_type;
    typedef real_type selection_type;
    typedef real_type dominance_type;

    typedef boost::random::uniform_01< real_type >  dist_type;
    typedef boost::random::bernoulli_distribution< real_type > neutral_dist_type;

    random_generator( URNG & rng ) :
        m_rng( &rng )
        , m_neutral()
    {}
    random_generator( URNG & rng, boost::property_tree::ptree & config ) : 
        m_rng( &rng )
        , m_neutral()
    {
        parseConfig( config );
    }

    random_generator( const self_type & other ) : 
        m_rng( other.m_rng )
        , m_neutral( other.m_neutral.param() )
    {}
    
    result_type operator()( ) {
        key_type k = m_uniform( *m_rng );

        selection_type s = DEFAULT_SELECTION;
        dominance_type d = DEFAULT_DOMINANCE;
        bool    is_neutral = m_neutral( *m_rng );
        if( !is_neutral ) {
            s = m_uniform( *m_rng );
            d = m_uniform( *m_rng );
        }

        return basic_allele(k, s, d, is_neutral );
    }

    virtual ~random_generator() {}
protected:
    void parseConfig( boost::property_tree::ptree & config ) {
        std::ostringstream oss;
        oss << CONFIG_BLOCK_K << "." << ALLELE_BLOCK_K << "." << NEUTRAL_P_K;
        if( config.get_child_optional( oss.str() ) == boost::none ) {
            config.put( oss.str(), m_neutral.p() );
        } else {
            double p = config.get< double >( oss.str(), m_neutral.p() );
            typename neutral_dist_type::param_type tmp( p );
            m_neutral.param(tmp);
        }
    }

    URNG *  m_rng;
    dist_type   m_uniform;
    neutral_dist_type m_neutral;
};

}   // namespace clotho {
}   // namespace utility {

#endif  // BASIC_ALLELE_GENERATOR_HPP_