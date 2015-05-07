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

#ifndef ASSORTATIVE_SELECTOR_HPP_
#define ASSORTATIVE_SELECTOR_HPP_

#include <boost/random/discrete_distribution.hpp>
#include <iterator>
#include <cmath>
#include <cassert>

template < class URNG >
struct assortative_selector {
    typedef URNG rng_type;
    typedef unsigned int result_type;

    typedef boost::random::discrete_distribution< result_type, double > dist_type;

    rng_type * m_rng;
    dist_type m_dist;

    unsigned int m_idx;

    assortative_selector( const assortative_selector< URNG > & other ) :
        m_rng( other.m_rng)
        , m_dist( other.m_dist )
        , m_idx( other.m_idx )
    {}

    template < class Iter >
    assortative_selector( rng_type & rng, Iter first, Iter last ) :
        m_rng( &rng )
        , m_dist()
        , m_idx( -1 )
    {
        
        unsigned int N = std::distance( first, last );

        N = ((N * ( N + 1 ) ) / 2);
        assert( N != -1 );

        std::vector< double > pairs( N, 0.0);
        std::vector< double >::iterator it = pairs.begin();
        Iter y = first;
        while( y != last ) {
            double fit = *y;
            ++y;
            Iter x = first;
            while( x != y ) {
                (*it) = fit * (*x);
                ++it;
                ++x;
                --N;
            }
        }

        assert( N == 0 );

        typename dist_type::param_type p(pairs.begin(), pairs.end());
        m_dist.param(p);
    }

    unsigned int operator()() {
        unsigned int ind = 0;
        if(m_idx == -1 ) {
            m_idx = m_dist( *m_rng );
            ind = (int) ((sqrt( 8 * m_idx + 1) - 1) / 2);
            m_idx = m_idx - (( ind + 1 ) * (ind ) / 2);
            // allows selfing
        } else {
            ind = m_idx;
            m_idx = -1;
        }
        return ind;
    }
};

#endif  // ASSORTATIVE_SELECTOR_HPP_
