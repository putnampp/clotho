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

#include "config_manager.hpp"

#include "simulate_engine.hpp"

#include "random_number_options.hpp"

int main( int argc, char ** argv ) {
    po::variables_map vm;

    int ret = config_manager::getInstance()->parse_commandline( argc, argv, vm );
    if(ret) return ret;

    typedef random_number_options::seed_type seed_type;

    seed_type seed = vm[ random_number_options::SEED_K ].as< seed_type >();

    simulate_engine eng( seed, 0.5, 0.5 );

    unsigned int i = 1;
    while( i-- ) {
        eng.simulate( 100 );

        std::cout << eng;
    }

    return 0;
}
