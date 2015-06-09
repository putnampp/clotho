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

#include <boost/program_options.hpp>
#include <boost/random/mersenne_twister.hpp>

#include "clotho/utility/timer.hpp"
#include "clotho/utility/popcount.hpp"

#include "square.h"
#include "square_stream.h"

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

template < class RNG, class GPUSquare >
void runGPUTest( RNG & rng, unsigned int list_len, unsigned int rounds, bool print ) {
    timer_type t;
    GPUSquare s( rng );
    t.stop();

    std::cout << "Init: " << t.elapsed().count() << std::endl;

    unsigned int i = 0;
    while( i < rounds ) {
        t.start();
        s(list_len);
        t.stop();

        std::cout << "Round " << (++i) <<  ": device random lapse: " << t.elapsed().count() << std::endl;
        if( print ) {
            std::cerr << "Round " << i << std::endl;
            std::cerr << s;
        }
    }
}

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
    std::cout << "host lapse:   " << t.elapsed().count() << std::endl;
}

template < class Engine >
void runCPUTest( Engine & eng, int_type N, unsigned int rounds, bool print ) {
    timer_type t;
    int_type * a = new int_type[ N ];
    t.stop();

    std::cout << "Init: " << t.elapsed().count() << std::endl;
    unsigned int i = 0;
    while( i < rounds ) {
        t.start();

        int_type * tmp = a;
        unsigned int j = N;
        while( j-- ) {
            (*tmp++) = popcount(eng());
        }

        t.stop();
        std::cout << "Round " << (++i) << ": host random lapse: " << t.elapsed().count() << std::endl;

        if( print ) {
            tmp = a;
            j = N;
            while( j-- ) {
                std::cerr << (N - j) << " -> " << (*tmp++) << std::endl;
            }
        }
    }

    delete [] a;
}

int main( int argc, char ** argv ) {

    po::variables_map vm;
    if( parse_commandline(argc, argv, vm ) ) {
        return -1;
    }

    rng_type rgen(vm[SEED_K].as<unsigned int>());
    unsigned int list_len = vm[LEN_K].as<unsigned int>();
    unsigned int rounds = vm[ROUNDS_K].as<unsigned int>();

    bool print_gpu = (vm.count(PRINT_GPU_K) > 0);
//    bool print_cpu = (vm.count(PRINT_CPU_K) > 0);


    switch( vm[TEST_GPU_K].as<unsigned int>() ) {
    case (unsigned int)-1:
        break;
    case 1:
        std::cout << "Testing CPU RNG -> GPGPU Squaring ... " << std::endl;
        runGPUTest<rng_type, SquareStream >( rgen, list_len, rounds, print_gpu );
        break;
    default:
        std::cout << "Testing GPGPU based RNG ..." << std::endl;
        runGPUTest<rng_type, Square >( rgen, list_len, rounds, false );
        break;
    }

    switch( vm[TEST_CPU_K].as<unsigned int>() ) {
    case (unsigned int)-1:
        break;
    default:
        std::cout << "Testing CPU based RNG ..." << std::endl;
        runCPUTest( rgen, list_len, rounds, false );
        break;
    }

    return 0;
}
