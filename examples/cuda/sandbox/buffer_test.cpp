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
//


/**
 * @section DESCRIPTION
 * 
 * This is a simple program that generates a stream of random numbers
 * using BufferedStream configurations
 *
 * The popcount of each random number is accumulated
 * The overall time for the main loop is recorded
 *
 * Four Tests are performed:
 * 1 - basic_buffer
 *   - Single-thread, refills buffer when buffer has been exhausted
 * 2 - double_buffered
 *   - Allocates 2*buffer-size space
 *   - Fills first half
 *   - Spawns worker thread to fill 'empty' half
 *   - Will block main thread until worker thread has completed by 
 *   -  joining the worker thread with main thread
 *   - (Scenario: repeatedly spawn worker thread)
 * 3 - double_buffered_mutex
 *   - Allocates 2*buffer-size space,
 *   - Fills first half
 *   - Spawns worker thread to fill 'empty' half
 *   - Uses mutex and conditional waits to fill 'empty' portion
 *   - Join of work thread occurs on object destruction
 *   - (Scenario: single worker pool; Used to evaluate thread creation overhead)
 * 4 - no_buffer
 *   - Performs the main loop without any threading
 *
 * The elapsed time for each test is recorded and a JSON formatted object is
 * printed to stdout
 */
#include <iostream>

#include <cuda.h>
#include <curand.h>

#include <boost/program_options.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "clotho/utility/timer.hpp"
#include "clotho/utility/popcount.hpp"
#include "clotho/utility/log_helper.hpp"

#include "buffered_stream.hpp"

namespace po=boost::program_options;

typedef clotho::utility::timer timer_type;

const std::string HELP_K = "help";

const std::string N_K = "length";
const std::string ROUND_K = "rounds";
const std::string BUFFER_SIZE_K = "buffer-size";
const std::string SEED_K = "seed";

class simple_generator : public igenerator {
public:
    typedef unsigned int result_t;

    simple_generator( unsigned long seed ) : igenerator( NULL, NULL ), m_gen(seed) {}

    void generate() {
        if( start == NULL || end == NULL ) return;

        result_t * t = static_cast< result_t *>( start );
        result_t * ee = static_cast< result_t *>( end );

        while( t < ee ) {
            *t++ = m_gen();
        }
    }

    virtual ~simple_generator() {}

protected:
    boost::random::mt19937 m_gen;
};

class cuda_generator : public igenerator {
public:
    typedef unsigned int result_t;

    cuda_generator( unsigned long seed ) : 
        igenerator( NULL, NULL )
        , m_dData( NULL )
        , m_size(0)
        , m_capacity(0)
        , m_seed_gen(seed)
        , m_gen()
        , m_status(true) {
        init();
    }

    
    void init( ) {
        curandStatus_t err =  curandCreateGenerator( &m_gen, CURAND_RNG_PSEUDO_MTGP32 );
        if( err != CURAND_STATUS_SUCCESS ) {
            std::cerr << "ERROR: Failed to create generator: " << err << std::endl;
            m_status = false;
            return;
        }
    }

    void generate() {
        if( cudaSetDevice(0) != cudaSuccess ) {
            std::cerr << "ERROR: Unable set CUDA device: " << std::endl;
            return;
        }

        if( !m_status ) { return; }

        resize();

        unsigned long long s = m_seed_gen();
        if( curandSetPseudoRandomGeneratorSeed( m_gen, s ) != CURAND_STATUS_SUCCESS ) {
            std::cerr << "ERROR: Unable to seed random number generator" << std::endl;
            m_status = false;
            return;
        }

        if( curandGenerate( m_gen, m_dData, m_size ) != CURAND_STATUS_SUCCESS ) {
            std::cerr << "ERROR: Unable to generate random numbers on device" << std::endl;
        }

        if( cudaMemcpy( start, m_dData, m_size * sizeof(result_t), cudaMemcpyDeviceToHost) != cudaSuccess) {
            std::cerr << "ERROR: Unable to copy data from device" << std::endl;
        }
    }
/*
    void generate_slow() {
        if( cudaSetDevice(0) != cudaSuccess ) {
            std::cerr << "ERROR: Unable set CUDA device: " << std::endl;
            return;
        }

        if( !m_status ) { return; }

        resize();

        curandGenerator_t g;
        curandStatus_t err =  curandCreateGenerator( &g, CURAND_RNG_PSEUDO_MTGP32 );
        if( err != CURAND_STATUS_SUCCESS ) {
            std::cerr << "ERROR: Failed to create generator: " << err << std::endl;
            m_status = false;
            return;
        }
     
        unsigned long long s = m_seed_gen();
        if( curandSetPseudoRandomGeneratorSeed( g, s ) != CURAND_STATUS_SUCCESS ) {
            std::cerr << "ERROR: Unable to seed random number generator" << std::endl;
            m_status = false;
            return;
        }

        if( curandGenerate( g, m_dData, m_size ) != CURAND_STATUS_SUCCESS ) {
            std::cerr << "ERROR: Unable to generate random numbers on device" << std::endl;
        }

        if( cudaMemcpy( start, m_dData, m_size * sizeof(result_t), cudaMemcpyDeviceToHost) != cudaSuccess) {
            std::cerr << "ERROR: Unable to copy data from device" << std::endl;
        }

        if( cudaDeviceSynchronize() != cudaSuccess ) {
            std::cerr << "ERROR: Sync error" << std::endl;
        }

        curandDestroyGenerator( g );
    }*/

    bool good() const { return m_status; }

    virtual ~cuda_generator() {
        if( m_status ) { curandDestroyGenerator( m_gen ); }

        if( m_dData ) { cudaFree( m_dData ); }
    }

protected:
    
    void resize( ) {
        size_t s = (result_t *)end - (result_t *)start;
//        std::cerr << "Size: " << s << " [" << start << ", " << end << "]" << std::endl;
        if( s > m_capacity ) {
            if( m_dData ) { cudaFree( m_dData ); }


            cudaDeviceSynchronize();
            cudaError_t err = cudaMalloc( (void **) &m_dData, s * sizeof( result_t ) );
            if( err != cudaSuccess ) {
                std::cerr << "ERROR: Unable to allocate device memory: " << std::endl;
                m_status = false;
            }
            m_capacity = s;
        }
        m_size = s;
    }

    result_t * m_dData;
    size_t     m_size, m_capacity;
    boost::random::mt19937 m_seed_gen;

    curandGenerator_t m_gen;
    bool m_status;
};

namespace clotho {
namespace utility {

/// BufferedStreams rely upon this helper to determine structure
template <>
struct result_type_of< simple_generator > {
    typedef unsigned int type;
};

template <>
struct result_type_of< cuda_generator > {
    typedef unsigned int type;
};

}   // namespace utility
}   // namespace clotho

template < class Func >
void runTimedFunc( Func & f, boost::property_tree::ptree & lapse ) {
    timer_type t;
    f();
    t.stop();
    clotho::utility::add_value_array( lapse, t );
}

template < class T >
struct summation {
    T val;

    summation( T v = 0) : val( v ) {}
    summation( const summation< T > & rhs ) : val(rhs.val) {}

    void operator()( T s ) { val += s; }
};

template < class T >
std::ostream & operator<<( std::ostream & out, const summation< T > & rhs ) {
    out << rhs.val;
    return out;
}

template < class Stream >
struct stream_summation {
    typedef typename clotho::utility::result_type_of< Stream >::type result_t;
    Stream * m_stream;
    unsigned int L;

    stream_summation( Stream & s, unsigned int l ) : m_stream( &s ), L(l) {}

    void operator()() {
        result_t res = 0;
        unsigned int n = L;
        while( n-- ) {
            res += (*m_stream)();   
        }
    }
};

template<>
struct stream_summation< void > {
    typedef unsigned int result_t;
    boost::random::mt19937 m_gen;
    unsigned int L;

    stream_summation( unsigned long seed, unsigned int l ) : m_gen(seed), L(l) {}

    void operator()() {
        result_t res = 0;
        unsigned int n = L;
        while( n-- ) {
            res += m_gen();
        }
    }
};

template < class Stream >
struct stream_popcount {
    typedef typename clotho::utility::result_type_of< Stream >::type result_t;
    Stream * m_stream;
    unsigned int L;

    boost::property_tree::ptree log;

    stream_popcount( Stream & s, unsigned int l ) : m_stream( &s ), L(l) {}

    void operator()() {
        result_t res = 0;
        unsigned int n = L;
        while( n-- ) {
            result_t tmp = (*m_stream)();
            res += popcount(tmp);

//            std::cerr << tmp << " = 0x" << std::hex << tmp << " -> " << std::dec << popcount(tmp) << std::endl;
            clotho::utility::add_value_array( log, tmp );
        }
    }
};

template<>
struct stream_popcount< void > {
    typedef unsigned int result_t;
    boost::random::mt19937 m_gen;
    unsigned int L;

    stream_popcount( unsigned long seed, unsigned int l ) : m_gen(seed), L(l) {}

    void operator()() {
        result_t res = 0;
        unsigned int n = L;
        while( n-- ) {
            res += popcount(m_gen());
        }
    }
};
void buildOptions( po::options_description & opts ) {
    po::options_description gen("General");

    gen.add_options()
    ((HELP_K + ",h").c_str(), "Print this")
    ;

    po::options_description params("Parameters");
    params.add_options()
    ((N_K + ",l").c_str(), po::value<unsigned int>()->default_value(100000), "Sequence Length")
    ((ROUND_K + ",r").c_str(), po::value<unsigned int>()->default_value(2), "Repeat tests")
    ((BUFFER_SIZE_K + ",b").c_str(), po::value<unsigned int>()->default_value(256), "BufferStream buffer size;  Size is number of elements to buffer")
    ((SEED_K + ",s" ).c_str(), po::value<unsigned long>()->default_value(1234), "Random Number Generator Seed")
    ;

    opts.add(gen).add(params);
}

int parse_cmd( int argc, char ** argv, po::variables_map & vm ) {
    po::options_description cmdline;
    buildOptions( cmdline );

    po::store( po::command_line_parser( argc, argv ).options(cmdline).run(), vm );

    int res = 0;
    if( vm.count( HELP_K ) ) {
        std::cout << cmdline << std::endl;
        res = 1;
    }

    return res;
}

int main( int argc, char ** argv ) {

    po::variables_map vm;

    if( parse_cmd( argc, argv, vm ) ) {
        return 1;
    }

#ifdef USE_CUDA_GENERATOR
    typedef cuda_generator                                              Generator;
#else
    typedef simple_generator                                            Generator;
#endif  // USE_CUDA_GENERATOR

    typedef BufferedStream< Generator, basic_buffer >                   SimpleBuffer;
    typedef BufferedStream< Generator, double_buffered >                DoubleBuffer;
    typedef BufferedStream< Generator, double_buffered_mutex >          DoubleBufferMutex;

    unsigned int rounds = vm[ROUND_K].as<unsigned int>();
    unsigned int N = vm[N_K].as<unsigned int>();
    unsigned long seed = vm[SEED_K].as<unsigned long>();
    unsigned int buffersize = vm[BUFFER_SIZE_K].as<unsigned int>();

    boost::random::mt19937 rng(seed);

    boost::property_tree::ptree lapses;

    std::cerr << "Basic Stream:" << std::endl;
    boost::property_tree::ptree basic_lapse, basic_total;
    for(unsigned int i = 0; i < rounds; ++i ) {
        timer_type t;
        Generator s(rng());
        SimpleBuffer simple_bs(s, buffersize);

        stream_popcount< SimpleBuffer > test( simple_bs, N );
        runTimedFunc( test, basic_lapse );
        t.stop();
        clotho::utility::add_value_array( basic_total, t );
    }
    lapses.add_child( "basic_stream.main", basic_lapse );
    lapses.add_child( "basic_stream.total", basic_total );
    lapses.put( "basic_stream.memory", buffersize * sizeof(SimpleBuffer::result_t));

    std::cerr << "Double Buffered Stream: " << std::endl;
    boost::property_tree::ptree double_lapse, double_total;
    for( unsigned int i = 0; i < rounds; ++i ) {
        timer_type t;
        Generator s2(rng());
        DoubleBuffer double_bs(s2, buffersize);

        stream_popcount< DoubleBuffer > test( double_bs, N );
        runTimedFunc(test, double_lapse);
        t.stop();
        clotho::utility::add_value_array( double_total, t );
    }
    lapses.add_child( "double_buffered_stream.main", double_lapse);
    lapses.add_child( "double_buffered_stream.total", double_total);
    lapses.put( "double_buffered_stream.memory", buffersize * sizeof(DoubleBuffer::result_t) * 2);

    std::cerr << "Double Buffered Mutex Stream:" << std::endl;
    boost::property_tree::ptree mutex_lapse, mutex_total, mutex_rand;
    for( unsigned int i = 0; i < rounds; ++i ) {
        timer_type t;
        Generator s2(rng());
        DoubleBufferMutex double_bs(s2, buffersize);

        stream_popcount< DoubleBufferMutex > test( double_bs, N );
        runTimedFunc(test, mutex_lapse);
        t.stop();
        clotho::utility::add_value_array( mutex_total, t );
        mutex_rand.push_back( std::make_pair("", test.log) );
    }
    lapses.add_child( "double_buffered_mutex_stream.main", mutex_lapse);
    lapses.add_child( "double_buffered_mutex_stream.total", mutex_total);
    lapses.put( "double_buffered_mutex_stream.memory", buffersize * sizeof(DoubleBufferMutex::result_t) * 2);
    lapses.add_child( "double_buffered_mutex_stream.random", mutex_rand );

    std::cerr << "Bufferless:" << std::endl;
    // assuming that initialization of this object is always negligible
    boost::property_tree::ptree no_lapse;
    for( unsigned int i = 0; i < rounds; ++i ) {
        stream_popcount< void > test( rng(), N );
        runTimedFunc(test, no_lapse);
    }
    lapses.add_child( "no_stream", no_lapse);

    lapses.put( "config.N", N );
    lapses.put( "config.seed", seed );
    lapses.put( "config.rounds", rounds );
    lapses.put( "config.buffersize", buffersize );

    boost::property_tree::write_json(std::cout, lapses );

    return 0;
}
