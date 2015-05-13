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
#include <boost/random/uniform_real_distribution.hpp>
#include <iterator>
#include <cmath>
#include <cassert>
#include <iostream>

template < class URNG >
struct assortative_selector {
    typedef URNG rng_type;
    typedef unsigned int result_type;

    typedef boost::random::discrete_distribution< result_type, double > dist_type;

    typedef std::vector< std::pair< double, unsigned int > > population_fitness_type;
    typedef typename population_fitness_type::iterator      fitness_iterator;

    static const unsigned int BIN_SIZE = 100;

    rng_type * m_rng;
    dist_type m_dist;
    boost::random::uniform_real_distribution< double > m_uniform;

    unsigned int m_next_bin;

    population_fitness_type m_pop_fit;

    struct simple_comp {
        bool operator()( const std::pair< double, unsigned int > & lhs, const std::pair< double, unsigned int > & rhs ) {
            return (lhs.first < rhs.first);
        }
    };

    assortative_selector( const assortative_selector< URNG > & other ) :
        m_rng( other.m_rng)
        , m_dist( other.m_dist )
        , m_uniform()
        , m_next_bin( other.m_next_bin )
        , m_pop_fit( other.m_pop_fit )
    {}

    template < class Iter >
    assortative_selector( rng_type & rng, Iter first, Iter last ) :
        m_rng( &rng )
        , m_dist()
        , m_uniform()
        , m_next_bin( -1 )
        , m_pop_fit()
    {
        unsigned int N = 0;
        while( first != last ) {
            m_pop_fit.push_back( std::make_pair( *first, N++ ) );
            ++first;
        }

        simple_comp comp;

        std::sort( m_pop_fit.begin(), m_pop_fit.end(), comp );

        std::vector< double > bin_mean;
        if( m_pop_fit.size() % BIN_SIZE ) {
            bin_mean.reserve( m_pop_fit.size() / BIN_SIZE + 1 );
        } else {
            bin_mean.reserve( m_pop_fit.size() / BIN_SIZE );
        }

        fitness_iterator pit = m_pop_fit.begin();
        double accum = 0.0;
        unsigned int b = BIN_SIZE;

        while( pit != m_pop_fit.end() ) {
            accum += pit->first;
            if( ! --b ) {
                bin_mean.push_back( accum / (double) BIN_SIZE );
                accum = 0.0;
                b = BIN_SIZE;
            }
            ++pit;
        }

        if( b != BIN_SIZE ) {
            bin_mean.push_back( accum / (BIN_SIZE - b) );
        }

        // N * (N + 1) / 2
        unsigned int bin_pairs = bin_mean.size() * (bin_mean.size() + 1)  / 2;
        assert( bin_pairs != -1 );

        std::vector< double > pairs( bin_pairs, 0.0);
        std::vector< double >::iterator it = pairs.begin();
        std::vector< double >::iterator y = bin_mean.begin();
        while( y != bin_mean.end() ) {
            double fit = *y;
            ++y;
            std::vector<double>::iterator x = bin_mean.begin();
            while( x != y ) {
                (*it) = fit * (*x);
                ++it;
                ++x;
            }
        }

        assert( it == pairs.end() );

        typename dist_type::param_type p(pairs.begin(), pairs.end());
        m_dist.param(p);
    }

    /// allows selfing
    unsigned int operator()() {
        double r = m_uniform( *m_rng );

        unsigned int ind = r * BIN_SIZE;

        if(m_next_bin == -1 ) {
            unsigned int bin_idx = m_dist( *m_rng );
            unsigned int bin_y = (int) ((sqrt( 8 * bin_idx + 1) - 1) / 2);
            m_next_bin = bin_idx - (( bin_y + 1 ) * (bin_y) / 2);

            ind = bin_y * BIN_SIZE + ind;
        } else {
            ind = m_next_bin * BIN_SIZE + ind;

            m_next_bin = -1;
        }

        if( ind >= m_pop_fit.size() ) {
            unsigned int last_bin_size = m_pop_fit.size() % BIN_SIZE;
            ind = m_pop_fit.size() - r * last_bin_size - 1;
        }

        ind = m_pop_fit[ind].second;
       
        return ind;
    }
};

#endif  // ASSORTATIVE_SELECTOR_HPP_
