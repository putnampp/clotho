#ifndef RECOMBINATION_GENERATOR_HPP_
#define RECOMBINATION_GENERATOR_HPP_

#include "clotho/utility/random_generator.hpp"
#include "clotho/recombination/recombination_def.hpp"

#include <boost/random/bernoulli_distribution.hpp>

extern const string BASE_SEQUENCE_BIAS_K;

namespace clotho {
namespace utility {

template < class URNG, class Sequence, class Classifier >
class random_generator< URNG, clotho::recombine::recombination< Sequence, Classifier > > {
public:
    typedef URNG                                                        rng_type;
    typedef clotho::recombine::recombination< Sequence, Classifier >    result_type;

    typedef clotho::utility::random_generator< URNG, Classifier >       classifier_generator_type;
    typedef typename classifier_generator_type::result_type             classifier_type;

    typedef double real_type;
    typedef boost::random::bernoulli_distribution< real_type >          dist_type;

    random_generator( rng_type & rng, boost::property_tree::ptree & config ) :
        m_rng( &rng )
        , m_cgen( rng, config )
        , m_dist( 0.5 )
    {
        parseConfig( config );
    }

    random_generator( rng_type & rng, classifier_generator_type & cgen, real_type p = 0.5 ) :
        m_rng( &rng )
        , m_cgen( cgen )
        , m_dist( p )
    {}

    result_type operator()() {
        classifier_type cfier = m_cgen();
        return result_type( cfier, m_dist( *m_rng ) );   
    }

protected:

    void parseConfig( boost::property_tree::ptree & config ) {
        std::ostringstream oss;
        oss << CONFIG_BLOCK_K << "." << REC_BLOCK_K << "." << BASE_SEQUENCE_BIAS_K;

        if( config.get_child_optional( oss.str() ) == boost::none ) {
            config.put( oss.str(), m_dist.p() );
        } else {
            double p = config.get< double >( oss.str(), 0.5 );
            typename dist_type::param_type tmp(p);
            m_dist.param( tmp );
        }
    }

    rng_type    * m_rng;
    classifier_generator_type m_cgen;
    dist_type   m_dist;
};

}   // namespace utility {
}   // namespace clotho {

#endif  // RECOMBINATION_GENERATOR_HPP_