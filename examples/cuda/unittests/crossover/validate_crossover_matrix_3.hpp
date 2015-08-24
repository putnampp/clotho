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
#ifndef VALIDATE_CROSSOVER_MATRIX_3_HPP_
#define VALIDATE_CROSSOVER_MATRIX_3_HPP_

#include "crossover_test.hpp"

#ifndef LOG_RANDOM_EVENTS
#define LOG_RANDOM_EVENTS
#endif  // LOG_RANDOM_EVENTS

#include "clotho/cuda/crossover/crossover_matrix_3.cuh"
#include "clotho/utility/log_helper.hpp"
#include "validate_event_sequence.hpp"

#include <algorithm>
#include <vector>
#include <unordered_map>

template < class HashVector, class EventMap >
void compute_prefix_sum( HashVector & h, EventMap & em ) {

    std::ostringstream oss;
    typename HashVector::value_type sum = 0;

    for( unsigned int i = 0; i < h.size(); ++i ) {
        typename HashVector::value_type p = h[i];
        h[i] = sum;
        sum += p;

        oss << "," << p;
    }

    h.push_back( sum );

    typename EventMap::iterator it = em.find( oss.str() );
    if( it == em.end() ) {
        em.insert( std::make_pair( oss.str(), 1 ) );
    } else {
        it->second++;
    }
}

template < class Iter, class HashVector, class EventVector >
unsigned int build_crossover( Iter first, Iter last, const HashVector & hash, const EventVector & evts ) {
    unsigned int res = 0, b = 1;
    while( first != last ) {
        int h = (int)( *first * (crossover< 3 >::THREADS_PER_BLOCK) );

        unsigned int lo = hash[ h++ ];
        unsigned int hi = hash[ h ];

        unsigned int c = lo;
        while( lo < hi ) {
            c += ( *first > evts[ lo++ ] );       
        }

        res |= (( c & 1) * b);

        ++first;
        b <<= 1;
    }
    return res;
}

bool validate( crossover_test< crossover< 3 > > &ct, boost::property_tree::ptree & err ) {

    typedef crossover< 3 > crossover_type;
    typedef crossover_test< crossover_type > test_type;

    typedef typename test_type::random_vector random_vector;
    
    typedef typename test_type::count_iterator count_iterator;
    typedef typename test_type::sequence_iterator sequence_iterator;
    typedef typename test_type::random_iterator random_iterator;
    typedef typename test_type::host_vector     host_vector;

    sequence_iterator s_it = ct.sequences.begin();
    count_iterator c_it = ct.event_list.begin();
    random_iterator e_it = ct.rand_pool.begin();

    unsigned int sequence_width = ct.allele_list.size() / crossover_type::ALLELE_PER_INT;


    typedef std::vector< typename crossover_type::event_count_type > hash_vector;
    typedef std::vector< typename crossover_type::real_type > event_vector;
    typedef std::vector< typename crossover_type::allele_type > allele_vector;
    typedef typename allele_vector::iterator    allele_iterator;

    bool valid = ((ct.allele_list.size() % crossover_type::ALLELE_PER_INT) == 0);
    if( !valid ) {
        err.put( "message", "Invalid number of alleles" );
        err.put( "observed", ct.allele_list.size() );
        err.put( "expected", ((ct.allele_list.size() / crossover_type::ALLELE_PER_INT) + 1) * crossover_type::ALLELE_PER_INT);
        return valid;
    }

    allele_vector alleles;
    alleles.reserve( ct.allele_list.size() );

    std::copy( ct.allele_list.begin(), ct.allele_list.end(), std::back_inserter( alleles ) );

    std::unordered_map< std::string, unsigned int > emap;

    int sidx = 0;
    while( valid ) {
        if( s_it == ct.sequences.end() ) {
            valid = ( c_it == ct.event_list.end() && e_it == ct.rand_pool.end() );

            if( !valid ) {
                err.put("message", "Reached end of sequences, but not end of event list or random pool");
            }
            break;
        } else if ( c_it == ct.event_list.end() || e_it == ct.rand_pool.end() ) {
            valid = false;
            err.put("message", "Reached end of either event list or random pool, but not end of sequences" );
            break;
        }

        hash_vector vHash;
        vHash.reserve( crossover_type::THREADS_PER_BLOCK + 1);

        event_vector vEvents;
        vEvents.reserve( crossover_type::THREADS_PER_BLOCK );

        std::copy( c_it, c_it + crossover_type::THREADS_PER_BLOCK, std::back_inserter( vHash ) );
        
        compute_prefix_sum( vHash, emap );
        
        std::copy( e_it, e_it + crossover_type::THREADS_PER_BLOCK, std::back_inserter( vEvents ) );

        valid = valid && std::is_sorted( vEvents.begin(), vEvents.begin() + vHash.back() );

        if( !valid ) {
            err.put("message", "Unordered event list" );
            err.put("sequence_index", sidx );

            boost::property_tree::ptree i, c, e;
            clotho::utility::add_value_array( i, c_it, c_it + crossover_type::THREADS_PER_BLOCK );
            clotho::utility::add_value_array( c, vHash.begin(), vHash.end() );
            clotho::utility::add_value_array( e, vEvents.begin(), vEvents.begin() + vHash.back() );

            err.add_child("events.input", i );
            err.add_child("events.prefix_sum", c );
            err.add_child("events.location", e );
            break;
        }

        c_it += crossover_type::THREADS_PER_BLOCK;
        e_it += crossover_type::THREADS_PER_BLOCK;
        allele_iterator a_it = alleles.begin();
        unsigned int soff = 0;
        while( a_it != alleles.end() ) {
            typename crossover_type::int_type   obs = (*s_it);
            typename crossover_type::int_type   exp = build_crossover( a_it, a_it + crossover_type::ALLELE_PER_INT, vHash, vEvents );

            valid = (obs == exp );
            if( !valid ) {
                err.put("message", "Unexpected crossover pattern" );
                err.put("sequence_index", sidx );
                err.put("sequence_offset", soff);

                std::ostringstream oss;
                oss << "0x" << std::hex << std::setfill('0') << std::setw( sizeof( crossover_type::int_type) * 2) << exp;

                err.put("expected", oss.str() );

                oss.str("");
                oss.clear();
                oss << "0x" << std::hex << std::setfill('0') << std::setw( sizeof( crossover_type::int_type) * 2) << obs;
                err.put("observed", oss.str() );
                
                boost::property_tree::ptree c, e, a;
                clotho::utility::add_value_array( c, vHash.begin(), vHash.end() );
                clotho::utility::add_value_array( e, vEvents.begin(), vEvents.end() );
                clotho::utility::add_value_array( a, a_it, a_it + crossover_type::ALLELE_PER_INT );

                err.add_child("events.prefix_sum", c );
                err.add_child("events.locations", e );
                err.add_child("alleles", a);
                break;
            }
            ++s_it;
            a_it += crossover_type::ALLELE_PER_INT;
            ++soff;
        }
        ++sidx;
    }

    valid = valid && validate_event_sequence_distribution( emap, err );

    return valid;
}

#endif  // VALIDATE_CROSSOVER_MATRIX_3_HPP_
