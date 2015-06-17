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
#include <boost/math/constants/constants.hpp>
#include <cmath>

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
    static const unsigned int BIN_COUNT = 100;

    rng_type * m_rng;
    dist_type m_dist;
    boost::random::uniform_real_distribution< double > m_uniform;

    unsigned int m_next_bin;

    population_fitness_type m_pop_fit;
    std::vector< size_t >   m_bin_sizes;

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
        , m_bin_sizes( other.m_bin_sizes )
    {}


    template < class Iter >
    assortative_selector( rng_type & rng, Iter first, Iter last ) :
        m_rng( &rng )
        , m_dist()
        , m_uniform()
        , m_next_bin( -1 )
        , m_pop_fit()
    {
        init_method2( first, last );
    }

    /// allows selfing
    unsigned int operator()() {
        unsigned int bstart = 0, bsize = 0;
        if( m_next_bin == -1 ) {
            unsigned int bin_idx = m_dist( *m_rng );
            unsigned int bin_y = (int) ((sqrt( 8 * bin_idx + 1) - 1) / 2);
            m_next_bin = bin_idx - (( bin_y + 1 ) * (bin_y) / 2);

            bstart = (( bin_y == 0 ) ? 0 : m_bin_sizes[ bin_y - 1] );
            bsize  = m_bin_sizes[ bin_y ] - bstart;
//            if( bsize == 0 ) {
//                std::cerr << bin_y << ", " << m_bin_sizes[bin_y] << std::endl;
//                std::cerr << "Bin Start: " << bstart << "; Bin Size: " << bsize << std::endl;
//                std::cerr << m_bin_sizes[bin_y] << std::endl;
//                std::cerr << m_dist.probabilities()[ bin_y ] << std::endl;
//                assert( false );
//            }
        } else {
            bstart = (( m_next_bin == 0 ) ? 0 : m_bin_sizes[ m_next_bin - 1] );
            bsize =  m_bin_sizes[m_next_bin] - bstart;
            m_next_bin = -1;
        }

        if( bsize == 0 ) {
            std::cerr << "Bin Start: " << bstart << "; Bin Size: " << bsize << std::endl;
            assert( false );
        }

        double r = m_uniform( *m_rng );
        return m_pop_fit[ bstart + r * bsize ].second;
    }

    virtual ~assortative_selector() {}
protected:

    template < class Iter >
    void init_method1( Iter first, Iter last ) {
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

        unsigned int bin_size_accum = 0;
        while( pit != m_pop_fit.end() ) {
            accum += pit->first;
            if( ! --b ) {
                bin_mean.push_back( accum / (double) BIN_SIZE );
                bin_size_accum += BIN_SIZE;
                m_bin_sizes.push_back( bin_size_accum );
                accum = 0.0;
                b = BIN_SIZE;
            }
            ++pit;
        }

        if( b != BIN_SIZE ) {
            bin_mean.push_back( accum / (BIN_SIZE - b) );
            bin_size_accum += (BIN_SIZE - b);
            m_bin_sizes.push_back( bin_size_accum );
        }
        assert(m_bin_sizes.back() == N );

        init_discrete_distribution( bin_mean );
    }

    template < class Iter >
    void init_method2( Iter first, Iter last ) {
//        std::cerr << "Init method 2 " << std::endl;
        unsigned int N = 0;
        double min_fit = 2.0, max_fit = 0.0;
        while( first != last ) {
            m_pop_fit.push_back( std::make_pair( *first, N++ ) );

            if( *first > max_fit ) {
                max_fit = *first;
            }

            if( *first < min_fit ) {
                min_fit = *first;
            }
            ++first;
        }

        simple_comp comp;
        std::sort( m_pop_fit.begin(), m_pop_fit.end(), comp );

        std::vector< double > bin_mean;
        bin_mean.reserve( BIN_COUNT );
        m_bin_sizes.reserve( BIN_COUNT );

//        double bin_size = (max_fit - min_fit) / (double) BIN_COUNT;
        if( max_fit == min_fit) { 
            // all individuals have the same fitness
            // therefore all bins have same average fitness
            // divide population into BIN_COUNT bins
                                                            // if N = 101 && BIN_COUNT = 100 then
            unsigned int ind_per_bin = (N / BIN_COUNT);     // ind_per_bin = 1
            unsigned int b = BIN_COUNT;                     // b = 100
            unsigned int n = N % BIN_COUNT;                 // n = 1
            unsigned int baccum = 0;

            while( n-- ) {                                  // single round
                // N - (BIN_COUNT - 1) * ind_per_bin 
                bin_mean.push_back( max_fit );
                baccum += (ind_per_bin + 1);
                m_bin_sizes.push_back( baccum ); // push_back(1 + 1)
                --b;                                        // b = 99
            }

            while( b-- ) {                                  // 99 rounds
                bin_mean.push_back( max_fit );              //
                baccum += (ind_per_bin);                    //
                m_bin_sizes.push_back( baccum );            // push_back( )
            }
                                                            // | m_bin_sizes | = | < 2, 3, 4, 5, ..., 100, 101 > | = 100
        } else {
            double bin_step = (1.0000001 * max_fit - min_fit) / (double)(BIN_COUNT);
            fitness_iterator pit = m_pop_fit.begin();
            double accum = 0.0, bin_max = min_fit + bin_step;

//            std::cerr << "Bin Range: <" << min_fit << ", " << max_fit << ", "  << bin_step << " >" << std::endl;

            unsigned int count = 0, bin_size_accum = 0;
            while( pit != m_pop_fit.end() ) {
                if( pit->first >= bin_max ) {
                
                    bin_mean.push_back( accum / (double) count );
                    bin_size_accum += count;
                    m_bin_sizes.push_back( bin_size_accum );

                    accum = 0.0;
                    count = 0;
                    bin_max += bin_step;

                    // advance to next containing bin
                    while( pit->first >= bin_max ) {
                        bin_mean.push_back( 0.0 );
                        m_bin_sizes.push_back( bin_size_accum );
                        bin_max += bin_step;
                    }
                }
                accum += pit->first;
                ++count;

                ++pit;
            }
            if( bin_size_accum < (double) N ) {
                bin_mean.push_back( accum / (double) count );
                bin_size_accum += count;
                assert( bin_size_accum == (double) N );

                m_bin_sizes.push_back( bin_size_accum );
            }
        }

//        std::cerr << "Last two bin sizes: " << m_bin_sizes[ m_bin_sizes.size() - 2] << ", " << m_bin_sizes.back() << std::endl;
        assert( m_bin_sizes.back() == N);
//        std::cerr << "bin count: " << bin_mean.size() << std::endl;
        assert( bin_mean.size() == BIN_COUNT );

        init_discrete_distribution( bin_mean );
    }

    void init_discrete_distribution( std::vector< double > & bin_mean ) {
        // N * (N + 1) / 2
//        std::cerr << "Bin Mean Size: " << bin_mean.size() << std::endl;
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
                (*it) = prob2( fit, (*x) );
                ++it;
                ++x;
            }
        }

        assert( it == pairs.end() );

        typename dist_type::param_type p(pairs.begin(), pairs.end());
        m_dist.param(p);
    }

    double prob1( double x, double y ) {
        return x * y;
    }

    double prob2( double x, double y, double mu = 0.0, double sigma = 1.0 ) {
        double res = x * y;
        if( res == 0.0 ) return 0.0;

        double coeff = (1.0/sigma) * boost::math::double_constants::one_div_root_two_pi;
        double var = sigma * sigma;

        res = ( (1.0 - res) - mu); // (1 - xy) - mu
        res *= res;
        res /= (2.0 * var);

        return coeff * exp(-res);
    }
};

#endif  // ASSORTATIVE_SELECTOR_HPP_
