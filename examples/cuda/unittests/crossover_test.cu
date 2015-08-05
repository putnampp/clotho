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
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <sstream>

#include <cuda.h>
#include <curand.h>
#include <thrust/device_vector.h>
#include <thrust/scan.h>

#include <boost/lexical_cast.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "clotho/cuda/curand_poisson_wrapper.hpp"
#include "clotho/cuda/curand_uniform_wrapper.hpp"

#include "clotho/cuda/crossover/crossover_matrix.cuh"

#include "clotho/utility/log_helper.hpp"

template < class CrossType >
struct crossover_test {
    typedef CrossType crossover_type;

    typedef thrust::device_vector< typename CrossType::real_type >          random_vector;
    typedef thrust::device_vector< typename CrossType::allele_type >        allele_vector;
    typedef thrust::device_vector< typename CrossType::event_count_type >   count_vector;
    typedef thrust::device_vector< typename CrossType::int_type >           sequence_vector;
    typedef std::vector< typename CrossType::real_type >                    host_vector; 

    typedef typename random_vector::iterator    random_iterator;
    typedef typename count_vector::iterator     count_iterator;
    typedef typename allele_vector::iterator    allele_iterator;
    typedef typename sequence_vector::iterator  sequence_iterator;
    typedef typename host_vector::iterator      host_iterator;

    random_vector       rand_pool;
    allele_vector       allele_list;
    count_vector        event_list;
    sequence_vector     sequences;

    template < class AlleleGenerator >
    void initialize( AlleleGenerator & aGen, size_t N ) {
        unsigned int tail = N % crossover_type::comp_cap_type::THREADS_PER_BLOCK;
        if( tail ) {
            N += (crossover_type::comp_cap_type::THREADS_PER_BLOCK - tail);
            std::cerr << "Warning: Increased Allele Size to - " << N << std::endl;
        }

        allele_list.resize(N);
        aGen( allele_list, N );
    }

    template < class CountGenerator, class EventGenerator >
    void simulate( CountGenerator & cGen, EventGenerator & eGen, size_t N ) {
        event_list.resize( N + 1 );
        cGen( event_list, N );
        thrust::exclusive_scan( event_list.begin(), event_list.end(), event_list.begin() );

        unsigned int nEvents = event_list.back();

        std::cout << "Events generated: " << nEvents << std::endl;

        rand_pool.resize( nEvents );

        eGen( rand_pool );

        unsigned int sequence_width = allele_list.size() / 32;

        sequences.resize( N * sequence_width );

        typename crossover_type::real_type * pool = rand_pool.data().get();
        typename crossover_type::allele_type * alleles = allele_list.data().get();
        typename crossover_type::event_count_type * events = event_list.data().get();
        typename crossover_type::int_type * seqs = sequences.data().get();

        crossover_type ct;

        ct( pool, alleles, events, seqs, N, allele_list.size(), sequence_width );
    }

    inline bool validate_empty_sequence( sequence_iterator first, sequence_iterator last ) {
        bool valid = true;
        while( valid && first != last ) {
            valid = (*first == 0);
            ++first;
        }
        return valid;
    }

    inline bool sort_randoms( host_vector & d, random_iterator first, random_iterator last ) {
        
        typename crossover_type::real_type sum = 0;
        while( first != last ) {
            sum -= log( *first );
            d.push_back( sum );
            ++first;
        }

        sum -= log( 0.1 );

        for( host_iterator it = d.begin(); it != d.end(); ++it ) {
            (*it) /= sum;
        }

        return true;
    }

    inline unsigned int determine_mask( host_vector & rec_events, allele_iterator first, allele_iterator last ) {
        unsigned int m = 0, bit = 1;
        while( first != last ) {
            host_iterator it = std::lower_bound( rec_events.begin(), rec_events.end(), *first );
            m |= (( rec_events.begin() - it ) & 1) * bit;
            bit <<= 1;
            ++first;
        }
        return m;
    }

    template < class Iter1, class Iter2 >
    bool validate_order( Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2 ) {
        const double epsilon = 0.000001;
        bool eq = true;
        while( eq ) {
            if( first1 == last1 ) {
                eq = (first2 == last2);
                break;
            } else if( first2 == last2 ) {
                eq = false;
                break;
            }

            eq = ( abs(*first1 - *first2) < epsilon );
            ++first1; ++first2;
        }
        return eq;
    }

    bool validate( boost::property_tree::ptree & err ) {
        std::cerr << "Validate" << std::endl;
        bool valid = true;

        const unsigned int ALLELE_PER_BLOCK = sizeof( typename crossover_type::int_type ) * 8;

        unsigned int sequence_width = allele_list.size() / ALLELE_PER_BLOCK;

        count_iterator c_it = event_list.begin();
        sequence_iterator s_it = sequences.begin();
        random_iterator r_it = rand_pool.begin();

        typename crossover_type::real_type * rands = rand_pool.data().get();

        random_vector ordered_events;

        unsigned int seq_idx = 0;
        while( valid && c_it + 1 != event_list.end() ) {
            std::cerr << "Validating sequence: " << seq_idx << std::endl;
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
                s_it += sequence_width;
                continue;
            }

            unsigned int nEvents = max_events - min_events;
            if( rand_pool.end() - r_it < nEvents ) {
                err.put("message", "Too few random numbers");
                err.put("sequence_index", seq_idx );
                err.put("expected", nEvents );
                err.put("observed", (rand_pool.end() - r_it));
                valid = false;
                break;
            }

            host_vector rec_events;
            valid = sort_randoms( rec_events, r_it, r_it + nEvents );
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
                clotho::utility::add_value_array( inp, r_it, r_it + nEvents);

                err.add_child( "expected", exp );
                err.add_child( "observed", obs );
                err.add_child( "input_sequence", inp );
                break;
            }

            r_it += nEvents;
            rands += nEvents;

            allele_iterator a_it = allele_list.begin();
            unsigned int nMasks = sequence_width;
            if( (sequences.end() - s_it) < sequence_width ) {
                err.put("message", "Unexpected sequence length" );
                err.put("sequence_index", seq_idx );
                err.put("expected", sequence_width );
                err.put("observed", (sequences.end() - s_it) );
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
};

int main(int argc, char ** argv ) {

    if( argc !=  6) {
        std::cerr << "<prog> <seed> <rho> <N> <samples> <alleles>" << std::endl;
        return 1;
    }

    boost::property_tree::ptree log;

    unsigned long long seed = boost::lexical_cast< unsigned long long >( argv[1] );
    crossover::real_type rho = boost::lexical_cast< crossover::real_type >( argv[2] );
    unsigned int N = boost::lexical_cast< unsigned int >(argv[3]);
    unsigned int samples = boost::lexical_cast< unsigned int >(argv[4] );
    unsigned int A = boost::lexical_cast< unsigned int >(argv[5]);

    log.put("parameters.seed", seed );
    log.put("parameters.rho", rho );
    log.put("parameters.N", N );
    log.put("parameters.samples", samples );
    log.put("parameters.alleles", A);
    
    curandGenerator_t dGen;
    if( curandCreateGenerator( &dGen, CURAND_RNG_PSEUDO_MTGP32 ) != CURAND_STATUS_SUCCESS ) {
        std::cerr << "Unable to create Generator" << std::endl;
        assert(false);
    }

    if( curandSetPseudoRandomGeneratorSeed( dGen, seed ) != CURAND_STATUS_SUCCESS ) {
        std::cerr << "Unable to set seed value" << std::endl;
        assert(false);
    }

    unsigned int s = samples;
    while( s-- ) {
        crossover_test< crossover > ct;
        clotho::cuda::fill_poisson< unsigned int, crossover::real_type> cGen( dGen, rho );
        clotho::cuda::fill_uniform< crossover::real_type > eGen(dGen);

        ct.initialize( eGen, A );
        ct.simulate( cGen, eGen, N );

        boost::property_tree::ptree err;
        if( !ct.validate(err) ) {
            std::ostringstream oss;
            oss << "sample." << (samples - s) << ".error";
            log.add_child(oss.str(), err );
        }
    }

    boost::property_tree::write_json(std::cout, log);

    return 0;
}
