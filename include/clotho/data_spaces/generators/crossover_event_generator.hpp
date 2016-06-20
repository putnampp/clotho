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
#ifndef CLOTHO_CROSSOVER_EVENT_GENERATOR_HPP_
#define CLOTHO_CROSSOVER_EVENT_GENERATOR_HPP_

#include <boost/property_tree/ptree.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <cstring>
#include <iostream>

#include "clotho/recombination/recombination_rate_parameter.hpp"
#include "clotho/recombination/sequence_bias_parameter.hpp"
#include "clotho/data_spaces/generators/position_distribution_helper.hpp"
#include "clotho/data_spaces/generators/crossover_event_distribution_helper.hpp"

#ifdef DEBUGGING
#define ASSERT_VALID_RANGE( x, min, max ) assert( min <= x && x < max );
#else
#define ASSERT_VALID_RANGE( x, min, max )
#endif

namespace clotho {
namespace genetics {

template < class RNG, class PositionType >
class crossover_event_generator {
public:
    typedef RNG             random_engine_type;
    typedef PositionType    position_type;
    typedef PositionType    event_type;
    typedef unsigned int    bin_type;

    typedef typename position_distribution_helper< PositionType >::type position_distribution_type;
    
    typedef crossover_event_distribution_helper< double > event_distribution_helper_type;

    typedef boost::random::bernoulli_distribution< double > sequence_bias_distribution_type;

    typedef typename event_distribution_helper_type::IntType    IntType;
    typedef typename event_distribution_helper_type::type       event_distribution_type;

    typedef position_type *                         position_vector;
    typedef event_type *                            event_vector;
    typedef bin_type *                              bin_vector;

    static const unsigned int BIN_MAX = 256;

    crossover_event_generator( random_engine_type * rng, boost::property_tree::ptree & config ) :
        m_rand( rng )
        , m_pos( NULL )
        , m_bins( NULL )
        , m_events( NULL )
        , m_events_size(0)
        , m_event_alloc(0)
        , m_base_seq(0)
        , m_bin_size(0)
        , m_alloc_size(0)
    {
        recombination_rate_parameter< double > rho( config);

        m_event_dist.param( typename event_distribution_type::param_type( rho.m_rho ) );

        sequence_bias_parameter< double > bias( config );
        m_seq_bias.param( typename sequence_bias_distribution_type::param_type( bias.m_bias ) );
    }

    void update( position_vector pos, size_t N ) {
        m_pos = pos;

        resizeBins( N );

        size_t M = N % 4;
        size_t i = 0;
        while( i < M ) {
            m_bins[ i ] = m_pos[i] * BIN_MAX;
            i += 1;
        }

        while( i < N ) {
            m_bins[ i ] = m_pos[i] * BIN_MAX;
            m_bins[ i + 1 ] = m_pos[ i + 1 ] * BIN_MAX;
            m_bins[ i + 2 ] = m_pos[ i + 2 ] * BIN_MAX;
            m_bins[ i + 3 ] = m_pos[ i + 3 ] * BIN_MAX;
            i += 4;
        }
    }

    // generate a new event distribution
    size_t generate() {
        // clear the counts
        IntType N = m_event_dist( *m_rand );    // generate the maximum number of events

        setBaseSequence( m_seq_bias( *m_rand ) );

        resize( N );

        size_t prev = 0;
        position_type accum = 0.0;

        size_t i = 0;
        position_type e = (position_type)(N);
        while( e > (position_type) 0 ) {
            position_type p = m_pos_dist( *m_rand );
            accum += log( p ) / e;

            p = (1.0 - exp( accum ));

            assert( i == 0 || m_events[ i - 1 ] < p);

            size_t cur = p * BIN_MAX;   // transform the current event to its specific bin
            while( prev < cur ) {
                m_counts[ prev++ ] = i;
            }

            m_events[ i++ ] = p;
            m_counts[ cur ] = i;
            prev = cur;
            e -= 1.0;
        }

        while( prev <= BIN_MAX ) {
            m_counts[ prev++ ] = i;
        }

        assert( i == N );
        return N;
    }

    // test whether the genetic position at the given index
    // 
    bool operator()( size_t index ) {
        ASSERT_VALID_RANGE( index, 0, m_bin_size )

        position_type p = m_pos[ index ];
        bin_type bin_index = m_bins[index];

        ASSERT_VALID_RANGE( bin_index, 0, BIN_MAX );

        size_t hi = m_counts[ bin_index ];
        size_t lo =  ((bin_index == 0) ? 0 : m_counts[ bin_index - 1]);

        while( lo < hi && m_events[lo] < p ){ ++lo; }

        return (m_base_seq ^ (lo & 1));    // % 2
    }

    void setBaseSequence( bool is_base ) {
        m_base_seq = ((is_base) ? 0 : 1);
    }

    size_t getBaseSequence() {
        return m_base_seq;
    }

    virtual ~crossover_event_generator() {
        if( m_events != NULL ) {
            delete [] m_events;
        }

        if( m_bins != NULL ) {
            delete [] m_bins;
        }
    }

protected:

    void resize( size_t e ) {
        if( e > m_event_alloc ) {
            if( m_events != NULL ) {
                delete [] m_events;
            }
            m_events = new event_type[ e ];

            m_event_alloc = e;
        }

        m_events_size = e;
    }

    void resizeBins( size_t N ) {
        if( N > m_alloc_size ) {
            if( m_bins != NULL ) {
                delete [] m_bins;
            }

            m_bins = new bin_type[ N ];
            m_alloc_size = N;
        }

        m_bin_size = N;
    }

    random_engine_type          * m_rand;
    event_distribution_type     m_event_dist;
    position_distribution_type  m_pos_dist;

    position_vector     m_pos;
    bin_vector          m_bins;

    event_vector        m_events;
    size_t              m_events_size, m_event_alloc;

    size_t              m_counts[ BIN_MAX + 1 ];

    sequence_bias_distribution_type m_seq_bias;
    size_t              m_base_seq;
    size_t              m_bin_size, m_alloc_size;
};

#undef ASSERT_VALID_RANGE

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_CROSSOVER_EVENT_GENERATOR_HPP_
