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
#ifndef VALIDATE_CROSSOVER_MATRIX_4_HPP_
#define VALIDATE_CROSSOVER_MATRIX_4_HPP_

#include "crossover_test.hpp"

#include <vector>

#include "validate_sequence.hpp"

#include "clotho/cuda/crossover/crossover_matrix_4.cuh"
#include "clotho/utility/log_helper.hpp"

template < class HashVector, class RandomVector, class ElementIterator >
unsigned int determine_mask( HashVector & hash, RandomVector & rands, ElementIterator first, ElementIterator last ) {
    typedef typename HashVector::value_type count_type;
    typedef typename RandomVector::value_type value_type;

    size_t bins = hash.size() - 1;
    double bin_size = (1.0 / (double)bins);
    unsigned int mask = 0, bit = 1;
    while( first != last ) {
        unsigned int bin_idx = (*first) * bins;
        double bin_start = ((double)bin_idx) * bin_size;
        assert( bin_idx < bins );

        count_type lo = hash[ bin_idx ];
        count_type hi = hash[ bin_idx + 1 ];

        value_type accum = 0;
        // find the first event
        while( lo < hi ) {
            accum += log( rands[ lo ] ) / ((value_type) (hi - lo));
            value_type event = bin_start + (1.0 - exp( accum )) * bin_size;
            if( event > *first ) {
                break;
            }
            ++lo;
        }

        mask |= ( (lo & 1) * bit );
        ++first;
        bit <<= 1;
    }
    return mask;
}

template < class DataVector, class Iter >
void build_bin_indexing( DataVector & b, Iter first, Iter last ) {
    while( first != last ) {
        b.push_back( *first * 1024 );
        ++first;
    }
}

template < class Iter1, class Iter2 >
bool validate_allele_indexing( Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2 ) {
    bool valid = true;

    while( valid ) {
        if( first1 == last1 ) {
            valid = (first2 == last2);
            break;
        } else if( first2 == last2 ) {
            valid = false;
            break;
        }

        valid = (*first1 == *first2);
        ++first1; ++first2;
    }

    return valid;
}

bool validate( crossover_test < crossover< 4 > > & ct, boost::property_tree::ptree & err ) {
    typedef crossover< 4 > crossover_type;
    typedef crossover_test< crossover_type > test_type;

    typedef typename test_type::random_vector random_vector;

    typedef typename test_type::count_iterator count_iterator;
    typedef typename test_type::sequence_iterator sequence_iterator;
    typedef typename test_type::random_iterator random_iterator;
    typedef typename test_type::host_vector     host_vector;

    bool valid = true;

    const unsigned int ALLELE_PER_BLOCK = sizeof( typename crossover_type::int_type ) * 8;

    unsigned int sequence_width = ct.allele_list.size() / ALLELE_PER_BLOCK;

    std::cerr << "Sequence Width: " << sequence_width << std::endl;
    std::cerr << "Event List Size: " << ct.event_list.size() << std::endl;

    count_iterator c_it = ct.event_list.begin();
    sequence_iterator s_it = ct.sequences.begin();
    random_iterator r_it = ct.rand_pool.begin();

    std::vector< int > expected_indexing;
    build_bin_indexing( expected_indexing, ct.allele_list.begin(), ct.allele_list.end() );

    thrust::device_vector< int > obs_indexing;
    obs_indexing.resize( ct.allele_list.size() );
    getAlleleIndex<<< dim3(1,1,1), dim3(32,32,1) >>>( ct.allele_list.data().get(), obs_indexing.data().get(), ct.allele_list.size() );

    valid = validate_allele_indexing( expected_indexing.begin(), expected_indexing.end(), obs_indexing.begin(), obs_indexing.end() );

    if( !valid ) {
        err.put("message", "Invalid allele bin indexing");

        boost::property_tree::ptree inp;
        clotho::utility::add_value_array( inp, ct.allele_list.begin(), ct.allele_list.end() );

        boost::property_tree::ptree exp;
        clotho::utility::add_value_array( exp, expected_indexing.begin(), expected_indexing.end() );

        boost::property_tree::ptree obs;
        clotho::utility::add_value_array( obs, obs_indexing.begin(), obs_indexing.end() );

        err.add_child( "allele.locations", inp );
        err.add_child( "bins.expected", exp );
        err.add_child( "bins.observed", obs );
    }

//    typename crossover_type::real_type * rands = ct.rand_pool.data().get();
//    random_vector ordered_events;

//    const unsigned int BIN_COUNT = 1024;
    const unsigned int BIN_COUNT = 32;
    typedef std::vector< typename crossover_type::event_count_type > event_hash_type;
    typedef typename event_hash_type::iterator hash_iterator;
    event_hash_type event_hash(BIN_COUNT + 1, 0);

    std::cerr << "Hash size: " << event_hash.size() << std::endl;

    unsigned int seq_idx = 0;
    while( valid && s_it != ct.sequences.end() ) {
//        std::cerr << seq_idx << std::endl;

        unsigned int i = 1;
        unsigned int nEvents = 0;
        while( i <= BIN_COUNT && c_it != ct.event_list.end() ) {
            nEvents += (*c_it);
            event_hash[i++] = nEvents;
            ++c_it;
        }

        if( nEvents == 0 ) {
            // no recombination events
            // corresponding sequence should be empty

            valid = validate_empty_sequence( s_it, s_it + sequence_width );
            if( !valid ) {
                err.put("message", "Sequence is not empty");
                err.put("sequence_index", seq_idx );
                break;
            }
            s_it += sequence_width;
            continue;
        }

        if( ct.rand_pool.end() - r_it < nEvents ) {
            err.put("message", "Too few random numbers");
            err.put("sequence_index", seq_idx );
            err.put("expected", nEvents );
            err.put("observed", (ct.rand_pool.end() - r_it));
            valid = false;
            break;
        }

        host_vector rec_events;
        rec_events.reserve(nEvents);

        int unused_rands = 100;
        while( nEvents-- ) {
            rec_events.push_back( *r_it );
            ++r_it;
            --unused_rands;
        }

        r_it += unused_rands;

//        valid = valid && (rec_events.size() == nEvents);
//        if( !valid ) {
//            err.put("message", "Ordering offset failed");
//            err.put("sequence_index", seq_idx );
//
//            boost::property_tree::ptree inp;
//            clotho::utility::add_value_array(inp, r_it, r_it + nEvents );
//
//            boost::property_tree::ptree obs;
//            clotho::utility::add_value_array(obs, rec_events.begin(), rec_events.end() );
//
//            err.add_child("input", inp );
//            err.add_child("observed", obs );
//            break;
//        }
//
//        ordered_events.resize( nEvents );
//
//        order_random<typename crossover_type::real_type, typename crossover_type::comp_cap_type ><<< dim3(1,1,1), dim3(32,32,1) >>>( rands, ordered_events.data().get(), nEvents );
//        cudaDeviceSynchronize();
//
//        valid = validate_order( rec_events.begin(), rec_events.end(), ordered_events.begin(), ordered_events.end() );
//        if( !valid ) {
//            err.put("message", "Unexpected event ordering" );
//            err.put("sequence_index", seq_idx );
//
//            // expected from host
//            boost::property_tree::ptree exp;
//            clotho::utility::add_value_array( exp, rec_events.begin(), rec_events.end() );
//
//            // observed from device
//            boost::property_tree::ptree obs;
//            clotho::utility::add_value_array( obs, ordered_events.begin(), ordered_events.end() );
//
//            boost::property_tree::ptree inp;
//            clotho::utility::add_value_array( inp, r_it, r_it + nEvents + 1);
//
//            err.add_child( "expected", exp );
//            err.add_child( "observed", obs );
//            err.add_child( "input_sequence", inp );
//            break;
//        }
//
//        r_it += nEvents + 1;
//        rands += nEvents + 1;

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
//            std::cerr << "bit block index: " << nMasks << std::endl;
            assert( a_it != ct.allele_list.end() );
            // observed from device
            typename crossover_type::int_type observed_mask = (*s_it);

            // expected from host
            typename crossover_type::int_type expected_mask = determine_mask(event_hash, rec_events, a_it, a_it + ALLELE_PER_BLOCK );

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

                boost::property_tree::ptree hsh;
                clotho::utility::add_value_array( hsh, event_hash.begin(), event_hash.end() );

                boost::property_tree::ptree rec;
                clotho::utility::add_value_array(rec, rec_events.begin(), rec_events.end());
                boost::property_tree::ptree all;
                clotho::utility::add_value_array(all, a_it, a_it + ALLELE_PER_BLOCK );

                err.add_child("recombination.hash", hsh);
                err.add_child("recombination.events", rec );
                err.put("recombination.size", rec_events.size() );
                err.add_child("allele.locations", all );
                
                break;
            }

            a_it += ALLELE_PER_BLOCK;
            ++s_it;
        }
        ++seq_idx;
    }

    return valid;
}
#endif  // VALIDATE_CROSSOVER_MATRIX_4_HPP_
