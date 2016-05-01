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
#include "qtlsim_config.hpp"

#include "engine_logger.hpp"

#include <iostream>
#include <algorithm>
#include <map>

#include <boost/random/mersenne_twister.hpp>

#include "qtlsim_engine.hpp"

typedef boost::random::mt19937          random_engine_type;
typedef Engine< random_engine_type >    simulate_engine_type;

int main( int argc, char ** argv ) {

    po::variables_map vm;
    int ret = config_manager_type::getInstance()->parse_commandline( argc, argv, vm );
    if( ret ) return ret;

    random_engine_type  rand;
    boost::property_tree::ptree config;

    simulate_engine_type engine( &rand, config );

    return 0;
}
