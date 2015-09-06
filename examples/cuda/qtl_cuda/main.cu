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

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "clotho/configuration_manager/config_manager.hpp"

//#include "clotho/cuda/allele_space/allele_space.hpp"
//#include "clotho/cuda/sequence_space/sequence_space.hpp"
//#include "clotho/cuda/population_space/population_space.hpp"
#include "clotho/cuda/data_spaces/data_space.hpp"

#include "options/configuration_option.hpp"

#include "qtl_cuda_simulate_engine.hpp"

typedef clotho::configuration_manager::config_manager       config_manager_type;

typedef float           real_type;
typedef unsigned int    int_type;

#ifdef USE_UNIT_ORDERING
typedef unit_ordered_tag< int_type > order_tag_type;
#else
typedef unordered_tag   order_tag_type;
#endif

typedef AlleleSpace< real_type, int_type, order_tag_type >  allele_space_type;
typedef SequenceSpace< int_type >                           sequence_space_type;

typedef PopulationSpace< sequence_space_type, allele_space_type > population_space_type;

typedef qtl_cuda_simulate_engine< population_space_type >   engine_type;

int main( int argc, char ** argv ) {

    po::variables_map vm;
    int ret = config_manager_type::getInstance()->parse_commandline( argc, argv, vm );
    if( ret ) return ret;

    boost::property_tree::ptree infile;
    bool print_config_only = true;
    if( vm.count( configuration_option::CONFIG_K ) ) {
        std::string p = vm[ configuration_option::CONFIG_K ].as< configuration_option::path_type >();

        boost::property_tree::read_json( p, infile );
        print_config_only = false;
    }

    boost::property_tree::ptree config = infile.get_child( "configuration", infile );
    
    engine_type sim_engine( config );

    if( !print_config_only ) {
        for( unsigned int i = 0; i < 100; ++i ) {
            sim_engine.simulate(10000);
        }
    }

    boost::property_tree::ptree ofile;
    ofile.add_child( "configuration", config );
    boost::property_tree::write_json( std::cout, ofile );

    return 0;
}
