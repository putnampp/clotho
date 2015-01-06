#ifndef INDIVIDUAL_SELECTOR_HPP_
#define INDIVIDUAL_SELECTOR_HPP_

#include <boost/random/discrete_distribution.hpp>

template < class URNG >
struct individual_selector {
    typedef URNG            rng_type;
    typedef unsigned int    result_type;

    typedef boost::random::discrete_distribution< result_type, double > dist_type;

    rng_type        * m_rng;
    dist_type       m_dist;


    individual_selector( const individual_selector< URNG > & other ) :
        m_rng( other.m_rng )
        , m_dist( other.m_dist.param() ) {
    }

    template < class Iter >
    individual_selector( rng_type & rng, Iter first, Iter last ) : m_rng( &rng ), m_dist( first, last ) {}

    unsigned int operator()() {
        return m_dist( *m_rng );
    }

    template < class Iter >
    void reset( Iter first, Iter last ) {
        typename dist_type::param_type p(first, last );
        m_dist.param( p );
    }
};

#endif  // INDIVIDUAL_SELECTOR_HPP_
