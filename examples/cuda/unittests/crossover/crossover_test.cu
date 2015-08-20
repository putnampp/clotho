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

//#include "clotho/cuda/crossover/crossover_matrix.cuh"
#include "crossover_test.hpp"
#include "crossover_validate.hpp"

#if USE_CROSSOVER_MATRIX == 2
#include "validate_crossover_matrix_2.hpp"
#elif USE_CROSSOVER_MATRIX == 3
#include "validate_crossover_matrix_3.hpp"
#elif USE_CROSSOVER_MATRIX == 4
#include "validate_crossover_matrix_4.hpp"
#else
#include "validate_crossover_matrix_5.hpp"
#endif

#include "clotho/utility/log_helper.hpp"
#include "clotho/utility/timer.hpp"
typedef clotho::utility::timer timer_type;

typedef crossover< USE_CROSSOVER_MATRIX > crossover_type;

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

//    typename base_type::real_type get_rate() {
//        return base_type::get_rate();
//    }
};

template < >
class event_generator_wrapper< crossover< 4 > > : public clotho::cuda::fill_poisson< unsigned int, crossover_type::real_type > {
public:
    typedef clotho::cuda::fill_poisson< unsigned int, crossover_type::real_type > base_type;
//    static const unsigned int BINS = 1024;
    static const unsigned int BINS = 32;

    boost::random::mt19937  m_rand;
    crossover_type::real_type    m_rate;

    event_generator_wrapper( curandGenerator_t g, crossover_type::real_type rate ) :
        base_type( g, rate / (double) BINS)
        , m_rand( clotho::utility::clock_type::now().time_since_epoch().count() )
        , m_rate( rate )
    {}

    unsigned int operator()( thrust::device_vector< unsigned int > & buf, size_t N ) {
        size_t eCount = BINS * N;
        buf.resize( eCount );

        ((base_type*)this)->operator()( buf, eCount );
        
//        boost::random::poisson_distribution< unsigned int, crossover_type::real_type > pDist( ((crossover_type::real_type )N ) * m_rate );
//        return pDist( m_rand );
        return N * 100;
    }
};

/*
template < >
class event_generator_wrapper< crossover< 3 > > {
public:

    event_generator_wrapper( curandGenerator_t g, crossover_type::real_type rate ) {}

    unsigned int operator()( thrust::device_vector< unsigned int > & buf, size_t N ) {
        return N;
    }
};*/

template < class CrossType >
struct pool_generator_wrapper {
    typedef clotho::cuda::fill_uniform< typename CrossType::real_type > type;
};

struct dummy_generator {
    dummy_generator( curandGenerator_t g ) {}

    template < class B >
    void operator()( B & buf, size_t N = 0 ) {}
};

template < >
struct pool_generator_wrapper< crossover< 3 > > {
    typedef dummy_generator type;
};


#if USE_CROSSOVER_MATRIX == 3
template <> template < class CountGenerator, class EventGenerator >
void crossover_test< crossover< 3 > >::simulate( CountGenerator & cGen, EventGenerator & eGen, size_t N ) {
    event_list.resize( N * 256 );
    rand_pool.resize( N * 256 );

    unsigned int sequence_width = allele_list.size() / crossover_type::ALLELE_PER_INT;

    sequences.resize( N * sequence_width );

    typename crossover_type::real_type * pool = rand_pool.data().get();
    typename crossover_type::allele_type * alleles = allele_list.data().get();
    typename crossover_type::event_count_type * events = event_list.data().get();
    typename crossover_type::int_type * seqs = sequences.data().get();

    ct( pool, alleles, events, seqs, N, allele_list.size(), sequence_width, cGen.get_rate() );

    cudaDeviceSynchronize();
}
#endif  // USE_CROSSOVER_MATRIX == 3

#if USE_CROSSOVER_MATRIX == 2
template<>    template < class AlleleGenerator >
void crossover_test< crossover< 2 > >::initialize( AlleleGenerator & aGen, size_t N ) {
    unsigned int tail = N % crossover_type::comp_cap_type::THREADS_PER_BLOCK;
    if( tail ) {
        N += (crossover_type::comp_cap_type::THREADS_PER_BLOCK - tail);
        std::cerr << "Warning: Increased Allele Size to - " << N << std::endl;
    }

    allele_list.resize(N);
    aGen( allele_list, N );

    ct.adjust_alleles( allele_list.data().get(), N );

    cudaDeviceSynchronize();
}

template <> template < class CountGenerator, class EventGenerator >
void crossover_test< crossover< 2 > >::simulate( CountGenerator & cGen, EventGenerator & eGen, size_t N ) {

    event_list.resize( N * THREAD_NUM );
    rand_pool.resize( N * 1024 );

    unsigned int sequence_width = allele_list.size() / crossover_type::ALLELE_PER_INT;

    sequences.resize( N * sequence_width );

    typename crossover_type::real_type * pool = rand_pool.data().get();
    typename crossover_type::allele_type * alleles = allele_list.data().get();
    typename crossover_type::event_count_type * events = event_list.data().get();
    typename crossover_type::int_type * seqs = sequences.data().get();

    ct( pool, alleles, events, seqs, N, allele_list.size(), sequence_width, cGen.get_rate() );

    cudaDeviceSynchronize();
}
#endif  // USE_CROSSOVER_MATRIX == 2

int main(int argc, char ** argv ) {

    if( argc !=  7) {
        std::cerr << "<prog> <seed> <rho> <N> <samples> <alleles> <skip_valid>" << std::endl;
        return 1;
    }

    boost::property_tree::ptree log;

    unsigned long long seed = boost::lexical_cast< unsigned long long >( argv[1] );
    crossover_type::real_type rho = boost::lexical_cast< crossover_type::real_type >( argv[2] );
    unsigned int N = boost::lexical_cast< unsigned int >(argv[3]);
    unsigned int samples = boost::lexical_cast< unsigned int >(argv[4] );
    unsigned int A = boost::lexical_cast< unsigned int >(argv[5]);
    bool skip_valid = boost::lexical_cast< bool >( argv[6] );

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

        event_generator_wrapper< crossover_type > cGen( dGen, rho );
        clotho::cuda::fill_uniform< typename crossover_type::real_type > aGen( dGen );   // allele generator
        typename pool_generator_wrapper< crossover_type >::type eGen(dGen);

        timer_type t;
        ct.initialize( aGen, A );
        t.stop();

        clotho::utility::add_value_array( init_perf_log, t );

        t.start();
        ct.simulate( cGen, eGen, N );
        t.stop();

        clotho::utility::add_value_array( sim_perf_log, t);

        boost::property_tree::ptree err;
        if( !skip_valid && !validate( ct, err ) ) {
            std::ostringstream oss;
            oss << "sample." << (samples - s);
            log.add_child(oss.str() + ".error", err );

            boost::property_tree::ptree _state;
            ct.get_state( _state );
            log.add_child( oss.str() + ".state", _state );
            break;
        }
    }

    log.add_child( "performance.initialization", init_perf_log );
    log.add_child( "performance.simulate", sim_perf_log );

    boost::property_tree::write_json(std::cout, log);

    return 0;
}
