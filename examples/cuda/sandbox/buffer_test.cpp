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
//    std::cout << "Lapse: " << t.elapsed().count() << std::endl;
    clotho::utility::add_value_array( lapse, t );
}

template < class Stream, class Func >
void runTest( Stream & s, Func f, unsigned int L ) {
    std::cout << "Before: " << f << std::endl;

    timer_type t;
    while( L-- ) {
        f( s() );
    }
    t.stop();

    std::cout << "Result: " << f << std::endl;
    std::cout << "Lapse: " << t.elapsed().count() << std::endl;
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

struct random_summation {
    typedef unsigned int result_t;
    boost::random::mt19937 m_gen;
    unsigned int L;

    random_summation( unsigned long seed, unsigned int l ) : m_gen(seed), L(l) {}

    void operator()() {
        result_t res = 0;
        unsigned int n = L;
        while( n-- ) {
            res += m_gen();
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

    unsigned int rounds = vm[ROUND_K].as<unsigned int>();
    unsigned int N = vm[N_K].as<unsigned int>();
    unsigned long seed = vm[SEED_K].as<unsigned long>();

    boost::random::mt19937 rng(seed);

    boost::property_tree::ptree basic_lapse;
    for(unsigned int i = 0; i < rounds; ++i ) {
        Generator s(rng(), 0);
        SimpleBuffer simple_bs(s, N/4);

        stream_summation< SimpleBuffer > test( simple_bs, N );
        runTimedFunc( test, basic_lapse );
    }

    boost::property_tree::ptree double_lapse;
    for( unsigned int i = 0; i < rounds; ++i ) {
        Generator s2(rng(), 0);
        DoubleBuffer double_bs(s2, N/5);

        stream_summation< DoubleBuffer > test( double_bs, N );
        runTimedFunc(test, double_lapse);
    }

    boost::property_tree::ptree no_lapse;
    for( unsigned int i = 0; i < rounds; ++i ) {
        random_summation test( rng(), N );
        runTimedFunc(test, no_lapse);
    }

    boost::property_tree::ptree lapses;
    lapses.add_child( "basic_stream", basic_lapse );
    lapses.add_child( "double_buffered_stream", double_lapse);
    lapses.add_child( "no_stream", no_lapse);

    boost::property_tree::write_json(std::cout, lapses );

    return 0;
}
