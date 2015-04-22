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
