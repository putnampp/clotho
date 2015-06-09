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

    simple_generator( unsigned long seed, unsigned int o = 0 ) : igenerator( NULL, NULL ), m_gen(seed) {}

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

namespace clotho {
namespace utility {

/// BufferedStreams rely upon this helper to determine structure
template <>
struct result_type_of< simple_generator > {
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

    stream_popcount( Stream & s, unsigned int l ) : m_stream( &s ), L(l) {}

    void operator()() {
        result_t res = 0;
        unsigned int n = L;
        while( n-- ) {
            res += popcount((*m_stream)());
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

    typedef simple_generator                                            Generator;
    typedef BufferedStream< Generator, basic_buffer >                   SimpleBuffer;
    typedef BufferedStream< Generator, double_buffered >                DoubleBuffer;
    typedef BufferedStream< Generator, double_buffered_mutex >          DoubleBufferMutex;

    unsigned int rounds = vm[ROUND_K].as<unsigned int>();
    unsigned int N = vm[N_K].as<unsigned int>();
    unsigned long seed = vm[SEED_K].as<unsigned long>();
    unsigned int buffersize = vm[BUFFER_SIZE_K].as<unsigned int>();

    boost::random::mt19937 rng(seed);

    boost::property_tree::ptree basic_lapse;
    for(unsigned int i = 0; i < rounds; ++i ) {
        Generator s(rng(), 0);
        SimpleBuffer simple_bs(s, buffersize);

        stream_popcount< SimpleBuffer > test( simple_bs, N );
        runTimedFunc( test, basic_lapse );
    }

    boost::property_tree::ptree double_lapse;
    for( unsigned int i = 0; i < rounds; ++i ) {
        Generator s2(rng(), 0);
        DoubleBuffer double_bs(s2, buffersize);

        stream_popcount< DoubleBuffer > test( double_bs, N );
        runTimedFunc(test, double_lapse);
    }

    boost::property_tree::ptree mutex_lapse;
    for( unsigned int i = 0; i < rounds; ++i ) {
        Generator s2(rng(), 0);
        DoubleBufferMutex double_bs(s2, buffersize);

        stream_popcount< DoubleBufferMutex > test( double_bs, N );
        runTimedFunc(test, mutex_lapse);
    }

    boost::property_tree::ptree no_lapse;
    for( unsigned int i = 0; i < rounds; ++i ) {
        stream_popcount< void > test( rng(), N );
        runTimedFunc(test, no_lapse);
    }

    boost::property_tree::ptree lapses;
    lapses.add_child( "basic_stream", basic_lapse );
    lapses.add_child( "double_buffered_stream", double_lapse);
    lapses.add_child( "double_buffered_mutex_stream", mutex_lapse);
    lapses.add_child( "no_stream", no_lapse);

    lapses.put( "config.N", N );
    lapses.put( "config.seed", seed );
    lapses.put( "config.rounds", rounds );
    lapses.put( "config.buffersize", buffersize );

    boost::property_tree::write_json(std::cout, lapses );

    return 0;
}
