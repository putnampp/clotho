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
#ifndef VALIDATE_CROSSOVER_MATRIX_5_HPP_
#define VALIDATE_CROSSOVER_MATRIX_5_HPP_

#include "crossover_test.hpp"

//#include <algorithm>
#include <sstream>

#include "validate_sequence.hpp"
#include "reorder_uniform.hpp"
#include "validate_order.hpp"
#include "determine_mask.hpp"

#include "clotho/cuda/crossover/crossover_matrix_5.cuh"
#include "clotho/utility/log_helper.hpp"

bool validate( crossover_test < crossover< 5 > > & ct, boost::property_tree::ptree & err ) {
    typedef crossover< 5 > crossover_type;
    typedef crossover_test< crossover_type > test_type;

    typedef typename test_type::random_vector random_vector;

    typedef typename test_type::count_iterator count_iterator;
    typedef typename test_type::sequence_iterator sequence_iterator;
    typedef typename test_type::random_iterator random_iterator;
    typedef typename test_type::host_vector     host_vector;

    bool valid = true;

    const unsigned int ALLELE_PER_BLOCK = sizeof( typename crossover_type::int_type ) * 8;

    unsigned int sequence_width = ct.allele_list.size() / ALLELE_PER_BLOCK;

    count_iterator c_it = ct.event_list.begin();
    sequence_iterator s_it = ct.sequences.begin();
    random_iterator r_it = ct.rand_pool.begin();

    typename crossover_type::real_type * rands = ct.rand_pool.data().get();

    random_vector ordered_events;

    unsigned int seq_idx = 0;
    while( valid && c_it + 1 != ct.event_list.end() ) {
        typename crossover_type::event_count_type min_events = *c_it;
        ++c_it;
        typename crossover_type::event_count_type max_events = *c_it;

        if( min_events >= max_events ) {
            // no recombination events
            // corresponding sequence should be empty

            valid = validate_empty_sequence( s_it, s_it + sequence_width );
            if( !valid ) {
                err.put("message", "Sequence is not empty");
                err.put("sequence_index", seq_idx );
                break;
            }
            ++r_it;
            ++rands;
            s_it += sequence_width;
            continue;
        }

        unsigned int nEvents = max_events - min_events;
        if( ct.rand_pool.end() - r_it < nEvents ) {
            err.put("message", "Too few random numbers");
            err.put("sequence_index", seq_idx );
            err.put("expected", nEvents );
            err.put("observed", (ct.rand_pool.end() - r_it));
            valid = false;
            break;
        }

        host_vector rec_events;
        
        //valid = sort_randoms( rec_events, r_it, r_it + nEvents + 1 );
        valid = reorder_uniform(r_it, r_it + nEvents + 1, std::back_inserter(rec_events) );

        valid = valid && (rec_events.size() == nEvents);
        if( !valid ) {
            err.put("message", "Ordering offset failed");
            err.put("sequence_index", seq_idx );

            boost::property_tree::ptree inp;
            clotho::utility::add_value_array(inp, r_it, r_it + nEvents );

            boost::property_tree::ptree obs;
            clotho::utility::add_value_array(obs, rec_events.begin(), rec_events.end() );

            err.add_child("input", inp );
            err.add_child("observed", obs );
            break;
        }

        ordered_events.resize( nEvents );

        order_random<typename crossover_type::real_type, typename crossover_type::comp_cap_type ><<< dim3(1,1,1), dim3(32,32,1) >>>( rands, ordered_events.data().get(), nEvents );
        cudaDeviceSynchronize();

        valid = validate_order( rec_events.begin(), rec_events.end(), ordered_events.begin(), ordered_events.end() );
        if( !valid ) {
            err.put("message", "Unexpected event ordering" );
            err.put("sequence_index", seq_idx );

            // expected from host
            boost::property_tree::ptree exp;
            clotho::utility::add_value_array( exp, rec_events.begin(), rec_events.end() );

            // observed from device
            boost::property_tree::ptree obs;
            clotho::utility::add_value_array( obs, ordered_events.begin(), ordered_events.end() );

            boost::property_tree::ptree inp;
            clotho::utility::add_value_array( inp, r_it, r_it + nEvents + 1);

            err.add_child( "expected", exp );
            err.add_child( "observed", obs );
            err.add_child( "input_sequence", inp );
            break;
        }

        r_it += nEvents + 1;
        rands += nEvents + 1;

        typename test_type::allele_iterator a_it = ct.allele_list.begin();
        unsigned int nMasks = sequence_width;
        if( (ct.sequences.end() - s_it) < sequence_width ) {
            err.put("message", "Unexpected sequence length" );
            err.put("sequence_index", seq_idx );
            err.put("expected", sequence_width );
            err.put("observed", (ct.sequences.end() - s_it) );
            valid = false;
            break;
        }

        while( nMasks-- ) {
            // observed from device
            typename crossover_type::int_type observed_mask = (*s_it);

            // expected from host
            typename crossover_type::int_type expected_mask = determine_mask( rec_events, a_it, a_it + ALLELE_PER_BLOCK );

            valid = (expected_mask == observed_mask);
            if( !valid ) {
                err.put("message", "Unexpected Mask" );
                err.put("sequence_index", seq_idx );
                err.put("block_index", (sequence_width - nMasks - 1));

                std::ostringstream oss;
                oss << "0x" << std::hex << std::setfill('0') << std::setw(8) << expected_mask;

                err.put("expected", oss.str() );

                oss.str("");
                oss.clear();
                oss << "0x" << std::hex << std::setfill('0') << std::setw(8) << observed_mask;

                err.put("observed", oss.str() );

                boost::property_tree::ptree rec;
                clotho::utility::add_value_array(rec, rec_events.begin(), rec_events.end());
                boost::property_tree::ptree all;
                clotho::utility::add_value_array(all, a_it, a_it + ALLELE_PER_BLOCK );
               
                err.add_child("recombination.events", rec );
                err.put("recombination.size", rec_events.size() );
                err.add_child("allele", all );
                break;
            }

            a_it += ALLELE_PER_BLOCK;
            ++s_it;
        }
        ++seq_idx;
    }

    return valid;
}
#endif  // VALIDATE_CROSSOVER_MATRIX_5_HPP_
