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
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/poisson_distribution.hpp>

#include "clotho/cuda/curand_poisson_wrapper.hpp"
#include "clotho/cuda/curand_uniform_wrapper.hpp"

//#include "crossover_test_5.h"
#include "clotho/cuda/crossover/crossover_matrix.cuh"
#include "crossover_test.hpp"
#include "crossover_validate.hpp"

#include "clotho/utility/log_helper.hpp"
#include "clotho/utility/timer.hpp"
typedef clotho::utility::timer timer_type;

typedef crossover< 5 > crossover_type;

template < class Method >
class event_generator_wrapper : public clotho::cuda::fill_poisson< unsigned int, crossover_type::real_type > {
public:
    typedef clotho::cuda::fill_poisson< unsigned int, crossover_type::real_type > base_type;

    event_generator_wrapper( curandGenerator_t g, crossover_type::real_type rate ) : base_type( g, rate ) {}

    unsigned int operator()( thrust::device_vector< unsigned int > & buf, size_t N ) {
        buf.resize( N + 1 );
        ((base_type*)this)->operator()( buf, N );

        thrust::exclusive_scan( buf.begin(), buf.end(), buf.begin() );

        return buf.back();
    }
};

template < unsigned int BinCount = 1024 >
struct event_hash_method {
    static const unsigned int BINS = BinCount;
};

template < unsigned int BC >
class event_generator_wrapper< event_hash_method< BC > > : public clotho::cuda::fill_poisson< unsigned int, crossover_type::real_type > {
public:
    typedef clotho::cuda::fill_poisson< unsigned int, crossover_type::real_type > base_type;
    typedef event_hash_method< BC > bin_type;

    boost::random::mt19937  m_rand;
    crossover_type::real_type    m_rate;

    event_generator_wrapper( curandGenerator_t g, crossover_type::real_type rate ) :
        base_type( g, rate / (double) bin_type::BINS)
        , m_rand( clotho::utility::clock_type::now().time_since_epoch().count() )
        , m_rate( rate )
    {}

    unsigned int operator()( thrust::device_vector< unsigned int > & buf, size_t N ) {
        size_t eCount = bin_type::BINS * N;
        buf.resize( eCount );

        ((base_type*)this)->operator()( buf, eCount );
        
        boost::random::poisson_distribution< unsigned int, crossover_type::real_type > pDist( ((crossover_type::real_type )N) * m_rate );

        return pDist( m_rand );
    }
};

/*
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
//        event_list.resize( N + 1 );
//        cGen( event_list, N );
//        thrust::exclusive_scan( event_list.begin(), event_list.end(), event_list.begin() );
//
//        unsigned int nEvents = event_list.back();

        unsigned int nEvents = cGen( event_list, N );
        std::cout << "Events generated: " << nEvents << std::endl;

        rand_pool.resize( nEvents + N );

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

//    inline bool validate_empty_sequence( sequence_iterator first, sequence_iterator last ) {
//        bool valid = true;
//        while( valid && first != last ) {
//            valid = (*first == 0);
//            ++first;
//        }
//        return valid;
//    }
//
//    inline bool sort_randoms( host_vector & d, random_iterator first, random_iterator last ) {
//        if( first == last ) return true;
//
//        typename crossover_type::real_type sum = 0;
//        while( true ) {
//            sum -= log( *first );
//            ++first;
//
//            if( first == last ) {
//                break;
//            }
//            d.push_back( sum );
//        }
//
//        for( host_iterator it = d.begin(); it != d.end(); ++it ) {
//            (*it) /= sum;
//        }
//
//        return true;
//    }

    inline unsigned int determine_mask( host_vector & rec_events, allele_iterator first, allele_iterator last ) {
        unsigned int m = 0, bit = 1;
        while( first != last ) {
            host_iterator it = std::upper_bound( rec_events.begin(), rec_events.end(), *first );
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

bool validate( crossover_test < crossover_matrix< 5 > > & ct, boost::property_tree::ptree & err ) {
    bool valid = true;

        const unsigned int ALLELE_PER_BLOCK = sizeof( typename crossover_type::int_type ) * 8;

        unsigned int sequence_width = allele_list.size() / ALLELE_PER_BLOCK;

        count_iterator c_it = ct.event_list.begin();
        sequence_iterator s_it = ct.sequences.begin();
        random_iterator r_it = ct.rand_pool.begin();

        typename crossover_type::real_type * rands = rand_pool.data().get();

        random_vector ordered_events;

        unsigned int seq_idx = 0;
        while( valid && c_it + 1 != event_list.end() ) {
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
            if( rand_pool.end() - r_it < nEvents ) {
                err.put("message", "Too few random numbers");
                err.put("sequence_index", seq_idx );
                err.put("expected", nEvents );
                err.put("observed", (rand_pool.end() - r_it));
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
};*/

int main(int argc, char ** argv ) {

    if( argc !=  6) {
        std::cerr << "<prog> <seed> <rho> <N> <samples> <alleles>" << std::endl;
        return 1;
    }

    boost::property_tree::ptree log;

    unsigned long long seed = boost::lexical_cast< unsigned long long >( argv[1] );
    crossover_type::real_type rho = boost::lexical_cast< crossover_type::real_type >( argv[2] );
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

    boost::property_tree::ptree init_perf_log, sim_perf_log;
    unsigned int s = samples;
    while( s-- ) {
        crossover_test< crossover_type > ct;
        //clotho::cuda::fill_poisson< unsigned int, crossover_type::real_type> cGen( dGen, rho );
        event_generator_wrapper< event_hash_method< 1024 > > cGen( dGen, rho );
        clotho::cuda::fill_uniform< crossover_type::real_type > eGen(dGen);

        timer_type t;
        ct.initialize( eGen, A );
        t.stop();

        clotho::utility::add_value_array( init_perf_log, t );

        t.start();
        ct.simulate( cGen, eGen, N );
        t.stop();

        clotho::utility::add_value_array( sim_perf_log, t);

        boost::property_tree::ptree err;
//        if( !ct.validate(err) ) {
        if( !validate( ct, err ) ) {
            std::ostringstream oss;
            oss << "sample." << (samples - s) << ".error";
            log.add_child(oss.str(), err );
            break;
        }
    }

    log.add_child( "performance.initialization", init_perf_log );
    log.add_child( "performance.simulate", sim_perf_log );

    boost::property_tree::write_json(std::cout, log);

    return 0;
}
