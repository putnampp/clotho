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
#ifndef CLOTHO_SMALL_CROSSOVER_EVENT_GENERATOR_HPP_
#define CLOTHO_SMALL_CROSSOVER_EVENT_GENERATOR_HPP_

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
class small_crossover_event_generator {
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

    typedef small_crossover_event_generator< RNG, PositionType > self_type;

    static const unsigned int BUFFER_MAX = 64;

    small_crossover_event_generator( random_engine_type * rng, double rho, double bias, position_vector pos, unsigned int N ) :
        m_rand( rng )
        , m_pos( pos )
        , m_pos_size( N )
        , m_events_size(0)
        , m_event_dist( rho )
        , m_seq_bias( bias )
    {}

    small_crossover_event_generator( random_engine_type * rng, boost::property_tree::ptree & config ) :
        m_rand( rng )
        , m_pos( NULL )
        , m_pos_size(0)
        , m_events_size(0)
        , m_base_seq(0)
    {
        recombination_rate_parameter< double > rho( config);

        m_event_dist.param( typename event_distribution_type::param_type( rho.m_rho ) );

        sequence_bias_parameter< double > bias( config );
        m_seq_bias.param( typename sequence_bias_distribution_type::param_type( bias.m_bias ) );
    }

    small_crossover_event_generator( const self_type & other ) :
        m_rand( other.m_rand )
        , m_pos( other.m_pos )
        , m_pos_size( other.m_pos_size )
        , m_events_size( other.m_events_size )
        , m_base_seq( other.m_base_seq )
    {   }

    void update( position_vector pos, size_t N ) {
        m_pos = pos;
        m_pos_size = N;
    }

    // generate a new event distribution
    size_t generate() {
        // clear the counts
        IntType N = m_event_dist( *m_rand );    // generate the maximum number of events
        resize( N );

        setBaseSequence( m_seq_bias( *m_rand ) );

        position_type accum = 0.0;

        size_t i = 0;
        position_type e = (position_type)(N);
        while( e > (position_type) 0 ) {
            position_type p = m_pos_dist( *m_rand );

            while( p == 0.0 ) { p = m_pos_dist( *m_rand ); }

            accum += log( p ) / e;

            p = (1.0 - exp( accum ));

            assert( i == 0 || m_events[ i - 1 ] < p);

            m_events[ i++ ] = p;
            e -= 1.0;
        }

        return N;
    }

    // test whether the genetic position at the given index
    // 
    inline bool operator()( const size_t index ) const {
#ifdef DEBUGGING
        if( index >= m_pos_size ) {
            std::cerr << "Crossover: " << index << " < " << m_pos_size << std::endl;
            assert(false);
        }
#endif

        position_type p = m_pos[ index ];
        return eval( p );
    }

    inline bool eval( position_type p ) const {

        bool res = false;
        if( m_events_size == 1 ) {
            res = ( m_events[0] < p );
        } else {
            size_t lo = 0;
            while( lo < m_events_size && m_events[lo] < p ){ ++lo; }
            res = ((lo & 1) ? true : false); // % 2
        }

        return res;
    }

    void setBaseSequence( bool is_base ) {
        m_base_seq = ((is_base) ? 0 : 1);
    }

    size_t getBaseSequence() const {
        return m_base_seq;
    }

    virtual ~small_crossover_event_generator() { }

protected:

    void resize( size_t e ) {
        assert( e < BUFFER_MAX );
        m_events_size = e;
    }

    random_engine_type          * m_rand;
    event_distribution_type     m_event_dist;
    position_distribution_type  m_pos_dist;

    position_vector     m_pos;
    size_t              m_pos_size;

    event_type          m_events[ BUFFER_MAX ];
    size_t              m_events_size;

    sequence_bias_distribution_type m_seq_bias;
    size_t              m_base_seq;
};

#undef ASSERT_VALID_RANGE

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_SMALL_CROSSOVER_EVENT_GENERATOR_HPP_
