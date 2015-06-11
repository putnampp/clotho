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
#include <stdio.h>
#include <cuda.h>
#include <iostream>
#include <cassert>
#include <sstream>

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include "clotho/utility/timer.hpp"
#include "clotho/utility/popcount.hpp"
#include "clotho/utility/log_helper.hpp"

#include "square.h"
#include "square_stream.h"
#include "square_host.h"

namespace po=boost::program_options;

typedef clotho::utility::timer  timer_type;
typedef boost::random::mt19937  rng_type;

typedef Square::int_type        int_type;

const std::string HELP_K = "help";

const std::string SEED_K = "seed";
const std::string LEN_K = "len";
const std::string ROUNDS_K = "rounds";
const std::string PRINT_GPU_K = "print-gpu";
const std::string PRINT_CPU_K = "print-cpu";

const std::string TEST_GPU_K = "test-gpu";
const std::string TEST_CPU_K = "test-cpu";

struct use_max {};

template < class RNG, class T = use_max >
struct param_gen {
    RNG * m_rng;
    unsigned int max;

    param_gen( RNG & r, unsigned int m ) : m_rng( &r ), max(m) {}

    unsigned int gen() {
        return max;
    }
};

struct use_uniform {};

template < class RNG >
struct param_gen< RNG, use_uniform > {
    RNG * m_rng;
    unsigned int max;

    boost::random::uniform_int_distribution <unsigned int> dist;

    param_gen( RNG & r, unsigned int m ) : m_rng( &r ), max( m ), dist(0, m) {}

    unsigned int gen() {
        return dist( *m_rng );
    }
};

/**
 * Given a Random Number Generator object, and Some object with the operator()( unsigned int ) defined
 *
 * A) Measure the time to construct an instance of object S from a pointer to RNG
 * B) For a number of rounds, measure the time it takes for S to operate given N (by default N is max)
 *
 */
template < class RNG, class S, class T = use_max >
struct test_func : public param_gen< RNG, T >  {
    typedef param_gen< RNG, T > base_type;
    boost::property_tree::ptree log;

/**
 * @param gen - random number generator
 * @param m - maximum value of N to be used in each round (by default N=max)
 */
    test_func( RNG & gen, unsigned int m = 1000000 ) : base_type( gen, m ) {}

    void operator()( unsigned int rounds ) {
        timer_type t;
        S obj( base_type::m_rng );
        t.stop();

        log.put("init", t.elapsed().count() );
        boost::property_tree::ptree l,d;

        while( rounds-- ) {
            unsigned int N = base_type::gen();

            t.start();
            obj( N );
            t.stop();
            clotho::utility::add_value_array(l, t);

            //std::ostringstream oss;
            //oss << obj;
            //clotho::utility::add_value_array(d, oss.str() );
        }

        log.add_child("rounds", l );
        log.add_child("data", d );
    }
};

template < class RNG, class T >
struct test_func< RNG, SquareHost, T > : public param_gen< RNG, T > {
    typedef param_gen< RNG, T > base_type;

    boost::property_tree::ptree log;

    test_func( RNG & gen, unsigned int m = 1000000 ) : base_type( gen, m ) {}

    void operator()( unsigned int rounds ) {
        timer_type t;
        SquareHost obj;
        t.stop();
        log.put("init", t.elapsed().count());

        boost::property_tree::ptree l,d;
        while( rounds-- ) {
            unsigned int N = base_type::gen();

            t.start();
            SquareHost::seed_type seed = (*base_type::m_rng)();
            obj( N, seed );
            t.stop();
            clotho::utility::add_value_array(l, t);

            assert( obj.good() );

            std::ostringstream oss;
            oss << obj;
            clotho::utility::add_value_array( d, oss.str() );
        }
        log.add_child("rounds", l );
        log.add_child("data", d );
    }
};

int parse_commandline( int argc, char ** argv, po::variables_map & vm ) {
    po::options_description gen( "General" );
    gen.add_options()
    ( (HELP_K + ",h").c_str(), "Print this" )
    ;

    po::options_description params( "Parameters" );
    params.add_options()
    ((SEED_K + ",s").c_str(), po::value< unsigned int >()->default_value(1234), "Host random number generator seed value" )
    ((LEN_K + ",l").c_str(), po::value< unsigned int >()->default_value( 2000 ), "List size" )
    ((ROUNDS_K + ",r").c_str(), po::value< unsigned int >()->default_value( 1 ), "Number of rounds to repeat")
    ((PRINT_GPU_K + ",p").c_str(), "Print GPGPU generated lists" )
    ((PRINT_CPU_K + ",P").c_str(), "Print CPU generated lists" )
    ((TEST_GPU_K + ",g").c_str(), po::value< unsigned int >()->default_value(0), "GPU test to run" )
    ((TEST_CPU_K + ",c").c_str(), po::value< unsigned int >()->default_value(0), "CPU test to run" )
    ;

    po::options_description cmdline;

    cmdline.add(gen).add(params);

    po::store( po::command_line_parser( argc, argv ).options(cmdline).run(), vm );

    int res = 0;
    if( vm.count( HELP_K ) ) {
        std::cout << cmdline << std::endl;
        res = 1;
    }

    return res;
}

void square_it( int_type N ) {

    timer_type t;
    int_type * a = new int_type[ N ];

    for(int_type i = 0; i < N; ++i ) {
        a[i] = i*i;
    }

    delete [] a;
    t.stop();
    std::cerr << "host lapse:   " << t.elapsed().count() << std::endl;
}

/**
 *
 * Resize buffer array if N is greater than current capacity
 *
 * Generate a random number
 * Calculate its popcount
 * Write value to array
 */
template < class RNG >
struct cpu_popcount {
    RNG * m_rng;
    unsigned int * m_buffer;
    unsigned int m_size, m_capacity;

    cpu_popcount( RNG & rng ) : m_rng(&rng), m_buffer(NULL), m_size(0), m_capacity(0) {}
    cpu_popcount( RNG * rng ) : m_rng(rng), m_buffer(NULL), m_size(0), m_capacity(0) {}

    void operator()( unsigned int N ) {
        if( N > m_capacity ) {
            if( m_buffer ) delete [] m_buffer;

            m_buffer = new unsigned int[ N ];
            m_capacity = N;
        }

        m_size = N;

        unsigned int * tmp = m_buffer;
        while( N-- ) {
            (*tmp++) = popcount( (*m_rng)() );       
        }
    }

    ~cpu_popcount() {
        if( m_buffer ) delete [] m_buffer;
    }
};

template < class R >
std::ostream & operator<<( std::ostream & out, const cpu_popcount<R> & rhs ) {
    if( rhs.m_size && rhs.m_buffer ) {
        unsigned int i = 0;
        out << rhs.m_buffer[ i++ ];       
        while( i < rhs.m_size ) {
            out << "," << rhs.m_buffer[i++];
        }
    }
    return out;
}

int main( int argc, char ** argv ) {

    po::variables_map vm;
    if( parse_commandline(argc, argv, vm ) ) {
        return -1;
    }

    rng_type rgen(vm[SEED_K].as<unsigned int>());
    unsigned int list_len = vm[LEN_K].as<unsigned int>();
    unsigned int rounds = vm[ROUNDS_K].as<unsigned int>();

//    bool print_gpu = (vm.count(PRINT_GPU_K) > 0);
//    bool print_cpu = (vm.count(PRINT_CPU_K) > 0);

    boost::property_tree::ptree log;

    switch( vm[TEST_GPU_K].as<unsigned int>() ) {
    case (unsigned int)-1:
        break;
    case 1: {
        std::cerr << "Testing CPU RNG -> GPGPU Squaring ... " << std::endl;
        test_func< rng_type, SquareStream > tf( rgen, list_len );
        tf( rounds );
        log.add_child( "cpu_stream", tf.log );
        break;
    }
    case 2: {
        std::cerr << "Testing GPGPU Host Squaring ..." << std::endl;
        test_func< rng_type, SquareHost > tf( rgen, list_len );
        tf( rounds );
        log.add_child( "square_host", tf.log );
        break;
    }
    default: {
        std::cerr << "Testing GPGPU based RNG ..." << std::endl;
        test_func< rng_type, Square > tf( rgen, list_len );
        tf( rounds );

        log.add_child( "gpu_only", tf.log );
        break;
    }
    }

    switch( vm[TEST_CPU_K].as<unsigned int>() ) {
    case (unsigned int)-1:
        break;
    default: {
        std::cerr << "Testing CPU based RNG ..." << std::endl;
        test_func< rng_type, cpu_popcount< rng_type > > tf( rgen, list_len );
        tf( rounds );

        log.add_child( "cpu_only", tf.log );
        break;
    }
    }

    boost::property_tree::write_json(std::cout, log );

    return 0;
}
