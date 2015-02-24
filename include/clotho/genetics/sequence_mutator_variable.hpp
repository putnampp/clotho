#ifndef SEQUENCE_MUTATOR_HPP_
#define SEQUENCE_MUTATOR_HPP_

#include "clotho/genetics/sequence_mutator_def.hpp"
#include "clotho/powerset/variable_subset.hpp"

template < class E, class B, class BM, class EK, class Generator >
class sequence_mutator< clotho::powersets::variable_subset< E, B, BM, EK >, Generator > {
public:
    typedef E           allele_type;
    typedef Generator   generator_type;

    typedef typename generator_type::result_type generator_result_type;

    typedef clotho::powersets::variable_subset< E, B, BM, EK > sequence_type;

    sequence_mutator( const generator_type & gen, unsigned int nEvents ) : m_gen(gen), m_events(nEvents) {}

    void operator()( std::shared_ptr< sequence_type > & seq ) {
        operator()( *seq );
    }

    void operator()( sequence_type & seq ) {
        for( unsigned int i = 0; i < m_events; ++i ) {
            generator_result_type a = m_gen( seq );
            seq.addElement( a );
        }
    }

    unsigned int event_count() const {
        return m_events;
    }

protected:
    generator_type   m_gen;
    unsigned int            m_events;
};

#include "clotho/utility/random_generator.hpp"
#include <boost/random/poisson_distribution.hpp>

namespace clotho {
namespace utility {

template < class URNG, class Sequence, class Generator >
class random_generator< URNG, sequence_mutator< Sequence, Generator > > {
public:
    typedef URNG                                        rng_type;
    typedef sequence_mutator< Sequence, Generator >     result_type;

    typedef Generator                                   mutation_generator_type;

    typedef boost::random::poisson_distribution< unsigned int, double > dist_type;

//    random_generator( rng_type & rng,  mutation_generator_type & mgen, double rate ) : m_rng(&rng), m_mgen( mgen ), m_dist(rate) {}
    random_generator( rng_type & rng, boost::property_tree::ptree & config ) :
        m_rng( &rng )
        , m_mgen( rng, config )
        , m_dist( DEFAULT_MUTATION_RATE ) {
        parseConfig( config );
    }

    result_type operator()() {
        unsigned int nEvents = m_dist(*m_rng);
        return result_type( m_mgen, nEvents );
    }

protected:

    void parseConfig( boost::property_tree::ptree & config ) {
        std::ostringstream oss;
        oss /*<< CONFIG_BLOCK_K << "."*/ << MUT_BLOCK_K << "." << RATE_PER_REGION_K;

        if( config.get_child_optional( oss.str() ) == boost::none ) {
            config.put( oss.str(), m_dist.mean() );
        } else {
            double mean = config.get< double >( oss.str(), DEFAULT_MUTATION_RATE );
            typename dist_type::param_type p( mean );
            m_dist.param( p );
        }
    }

    rng_type    * m_rng;
    mutation_generator_type m_mgen;
    dist_type   m_dist;
};

}   // namespace utility {
}   // namespace clotho {

#endif  // SEQUENCE_MUTATOR_HPP_
