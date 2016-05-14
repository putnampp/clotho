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
#include <cstring>
#include <iostream>

#include "clotho/recombination/recombination_rate_parameter.hpp"
#include "clotho/data_spaces/generators/position_distribution_helper.hpp"
#include "clotho/data_spaces/generators/crossover_event_distribution_helper.hpp"

namespace clotho {
namespace genetics {

template < class RNG, class PositionType >
class crossover_event_generator {
public:
    typedef RNG             random_engine_type;
    typedef PositionType    position_type;
    typedef PositionType    event_type;

    typedef typename position_distribution_helper< PositionType >::type position_distribution_type;
    
    typedef crossover_event_distribution_helper< double > event_distribution_helper_type;

    typedef typename event_distribution_helper_type::IntType    IntType;
    typedef typename event_distribution_helper_type::type       event_distribution_type;

    typedef std::vector< position_type >                position_vector;
    typedef typename position_vector::iterator          position_iterator;
    typedef typename position_vector::const_iterator    const_position_iterator;

    typedef event_type *                                event_vector;

    typedef std::vector< size_t >                       bin_vector;

    static const unsigned int BIN_MAX = 256;

    crossover_event_generator( random_engine_type * rng, boost::property_tree::ptree & config ) :
        m_rand( rng )
        , m_events( NULL )
        , m_events_size(0)
    {
        recombination_rate_parameter< double > rho( config);

        m_event_dist.param( typename event_distribution_type::param_type( rho.m_rho ) );
    }

    // copy the genetic position of all alleles into a local vector
    void update( position_iterator first, position_iterator last ) {
        m_pos.clear();
        m_bins.clear();
        while( first != last ) {
            position_type p = *first++;

            m_pos.push_back(p);

            assert( p < 1.0 );

            assert( p * BIN_MAX < BIN_MAX );
            m_bins.push_back( p * BIN_MAX );
        }
    }

    // generate a new event distribution
    void generate() {
        // clear the counts
        IntType e = m_event_dist( *m_rand );    // generate the maximum number of events

        resize( e );

        size_t prev = 0;
        position_type accum = 0.0;

        size_t i = 0;
        while( e ) {
            position_type p = m_pos_dist( *m_rand );
            accum += log( p ) / (position_type) e;

            p = (1.0 - exp( accum ));

            assert( i == 0 || m_events[ i - 1 ] < p);

            size_t cur = p * BIN_MAX;   // transform the current event to its specific bin
            while( prev < cur ) {
                m_counts[ prev++ ] = i;
            }

            m_events[ i++ ] = p;
            m_counts[ cur ] = i;
            prev = cur;
            --e;
        }

        while( prev <= BIN_MAX ) {
            m_counts[ prev++ ] = i;
        }

    }

    // test whether the genetic position at the given index
    // 
    bool operator()( size_t index ) {

        assert( 0 <= index && index < m_pos.size() );

        position_type p = m_pos[ index ];
        size_t bin_index = m_bins[index];

        assert( bin_index < BIN_MAX );

        size_t hi = m_counts[ bin_index ];
        size_t lo =  ((bin_index == 0) ? 0 : m_counts[ bin_index - 1]);

        while( lo < hi && m_events[lo] < p ){ ++lo; }

        return (lo & 1);    // % 2
    }

    virtual ~crossover_event_generator() {
        if( m_events != NULL ) {
            delete [] m_events;
        }
    }

protected:

    void resize( size_t e ) {
        if( e > m_events_size ) {
            if( m_events != NULL ) {
                delete [] m_events;
            }
            m_events = new event_type[ e ];

            m_events_size = e;
        }
    }

    random_engine_type          * m_rand;
    event_distribution_type     m_event_dist;
    position_distribution_type  m_pos_dist;

    position_vector     m_pos;
    bin_vector          m_bins;

    event_vector        m_events;
    size_t              m_events_size;

    size_t              m_counts[ BIN_MAX + 1 ];
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_CROSSOVER_EVENT_GENERATOR_HPP_
