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
#ifndef VALIDATE_CROSSOVER_MATRIX_2_HPP_
#define VALIDATE_CROSSOVER_MATRIX_2_HPP_

#include "crossover_test.hpp"

#ifndef LOG_RANDOM_EVENTS
#define LOG_RANDOM_EVENTS
#endif  // LOG_RANDOM_EVENTS

#include "clotho/cuda/crossover/crossover_matrix_2.cuh"

#include <cassert>
#include <vector>

#include "clotho/utility/log_helper.hpp"

template < class Iter >
bool validate_allele_ordering( Iter first, Iter last, boost::property_tree::ptree & err ) {
    bool valid = true;

    unsigned int idx = 0;
    while( valid && first != last ) {
        unsigned int lane_id = (unsigned int) (*first * 32.0);

        valid = ((idx & 31) == lane_id);
        if( ! valid ) {
            err.put( "message", "Unexpected lane id for allele" );

            err.put( "allele.index", idx );
            err.put( "allele.location", *first );
            err.put( "lane_id.expected", (idx & 31) );
            err.put( "lane_id.observed", lane_id );
            break;
        }
        ++idx;
        ++first;
    }

    return valid;
}

template < class Iter >
std::vector< unsigned int > build_event_range( Iter first, Iter last ) {
    unsigned int max = 0;
    unsigned int sum = 0;
    std::vector< unsigned int > res;
    while( first != last ) {
        if( *first > max ) {
            max = *first;
        }
        res.push_back( sum );
        sum += *first;
        ++first;
    }

    unsigned int N = res.size();
    for( unsigned int i = 1; i < N; ++i ) {
        res.push_back( res[i] );
    }
    res.push_back( sum );

    assert( res.size() == 2 * N );

    for( unsigned int i = 0; i < N; ++i ) {
        res.push_back( max );
    }

    assert( res.size() == 3 * N );

    return res;
}

template < class EvtVector >
bool validate_event_range( EvtVector & counts, boost::property_tree::ptree & err ) {
    if( (counts.size() & (THREAD_NUM - 1))  != 0 ) {
        err.put("message", "Unexpected size of event count list");
        err.put("list_size.observed", counts.size() );
        err.put("list_size.expected", (counts.size() & (~(THREAD_NUM - 1))) + THREAD_NUM);
        return false;
    }

    bool valid = true;
    EvtVector obs;
    obs.resize( 3 * 32 );

    typename EvtVector::value_type * c_ptr = counts.data().get();

    typedef typename EvtVector::iterator device_iterator;
    device_iterator c_it = counts.begin();

    typedef std::vector< typename EvtVector::value_type > host_vector;
    typedef typename host_vector::iterator host_iterator;

    unsigned int s_idx = 0;
    while( valid && c_it != counts.end() ) {
        get_warp_event_range<<< 1, 32 >>>( c_ptr, obs.data().get() );
        host_vector exp = build_event_range( c_it, c_it + 32 );
        cudaDeviceSynchronize();

        host_iterator h_it = exp.begin();
        device_iterator d_it = obs.begin();
        unsigned int e_idx = 0;
        while( true ) {
            if( h_it == exp.end() ) {
                valid = (d_it == obs.end() );

                if( !valid ) {
                    err.put( "message", "Unexpected early end to expected event range data");
                    err.put( "sequence.index", s_idx );
                    boost::property_tree::ptree o, e;
                    clotho::utility::add_value_array( o, obs.begin(), obs.end() );
                    clotho::utility::add_value_array( e, exp.begin(), exp.end() );

                    err.add_child( "event_range.observed", o );
                    err.add_child( "event_range.expected", e );
                }
                break;
            } else if( d_it == obs.end() ) {
                valid = false;

                err.put("message", "Unexpected early end to observed event range data" );
                err.put("sequence.index", s_idx );

                boost::property_tree::ptree o, e;
                clotho::utility::add_value_array( o, obs.begin(), obs.end() );
                clotho::utility::add_value_array( e, exp.begin(), exp.end() );

                err.add_child( "event_range.observed", o );
                err.add_child( "event_range.expected", e );
                break;
            }

            valid = ( *d_it == *h_it );

            if( !valid ) {
                if( e_idx < 32 ) {
                    err.put("message", "Mismatch in event low range" );
                } else if( e_idx < 64 ) {
                    err.put( "message", "Mismatch in event high range" );
                } else {
                    err.put( "message", "Mismatch in event maximum" );
                }

                err.put("sequence.index", s_idx );
                err.put("range.err_index", e_idx );

                
                boost::property_tree::ptree i, o, e;
                clotho::utility::add_value_array( i, c_it, c_it + 32 );
                clotho::utility::add_value_array( o, obs.begin(), obs.end() );
                clotho::utility::add_value_array( e, exp.begin(), exp.end() );

                err.add_child( "range.input", i );
                err.add_child( "range.expected", e );
                err.add_child( "range.observed", o );
                break;
            }

            ++d_it;
            ++h_it;
            ++e_idx;
        }

        c_ptr += THREAD_NUM;  
        c_it += THREAD_NUM;
        ++s_idx;
    }

    return valid;
}

bool validate( crossover_test< crossover< 2 > > & ct, boost::property_tree::ptree & err ) {
    typedef crossover< 2 > crossover_type;
    typedef crossover_test< crossover_type > test_type;

    typedef typename test_type::random_vector random_vector;
    typedef typename test_type::count_vector count_vector; 
   
    typedef typename test_type::count_iterator count_iterator;
    typedef typename test_type::sequence_iterator sequence_iterator;
    typedef typename test_type::random_iterator random_iterator;
    typedef typename test_type::host_vector     host_vector;

    // validate that alleles are ordered by lane
    bool valid = validate_allele_ordering( ct.allele_list.begin(), ct.allele_list.end(), err );

    valid = valid && validate_event_range( ct.event_list, err );

    return valid;
}

#endif  // VALIDATE_CROSSOVER_MATRIX_2_HPP_
