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
#include "mutation_rate_option.hpp"
#include "recombination_rate_option.hpp"
#include "population_size_option.hpp"

#include "clotho/utility/timer.hpp"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

int main( int argc, char ** argv ) {
    po::variables_map vm;

    int ret = config_manager::getInstance()->parse_commandline( argc, argv, vm );
    if(ret) return ret;

    typedef random_number_options::seed_type seed_type;
    seed_type seed = vm[ random_number_options::SEED_K ].as< seed_type >();

    typedef mutation_rate_option::mutation_rate_type mutation_rate_type;
    mutation_rate_type mu = vm[ mutation_rate_option::RATE_K ].as< mutation_rate_type >();

    typedef recombination_rate_option::recombination_rate_type recombination_rate_type;
    recombination_rate_type rho = vm[ recombination_rate_option::RATE_K ].as< recombination_rate_type >();

    typedef population_size_option::population_size_type population_size_type;
    population_size_type pop_size = vm[ population_size_option::SIZE_K ].as< population_size_type >();

    simulate_engine eng( seed, mu, rho, pop_size );

    boost::property_tree::ptree log;

    unsigned int i = 3;
    while( i-- ) {

        clotho::utility::timer t;
        eng.simulate( pop_size );
        t.stop();

//        boost::property_tree::ptree n;
//        eng.record_state( n );
//
//        std::ostringstream oss;
//        oss << "log" << i;
//        log.add_child( oss.str(), n );

        std::cerr << "Simulate Lapse: " << t << std::endl;
    }
/*
    test_log tl;
    test_log2 ts( &tl );

    ts.log_state( log, eng );*/

//    boost::property_tree::write_json( std::cout, log );

    return 0;
}
