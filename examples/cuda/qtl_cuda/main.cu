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
#include "clotho/genetics/population_growth_toolkit.hpp"

#include "simulation_log.hpp"

typedef clotho::configuration_manager::config_manager       config_manager_type;

#ifdef USE_DOUBLE_REAL
typedef double          real_type;
#else
typedef float           real_type;
#endif
typedef unsigned int    int_type;

#ifdef USE_UNIT_ORDERING
typedef unit_ordered_tag< int_type > order_tag_type;
#else
typedef unordered_tag   order_tag_type;
#endif

typedef AlleleSpace< real_type, int_type, order_tag_type >  allele_space_type;
typedef SequenceSpace< int_type >                           sequence_space_type;

//typedef PopulationSpace< sequence_space_type, allele_space_type > population_space_type;
typedef PopulationSpace< real_type, int_type, order_tag_type > population_space_type;

typedef qtl_cuda_simulate_engine< population_space_type >   engine_type;

static const std::string GEN_K = "generations";

typedef std::shared_ptr< ipopulation_growth_generator >                     population_growth_generator_type;
typedef std::shared_ptr< ipopulation_growth >                               population_growth_type;

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

    unsigned int nGens = 100;

    if( config.get_child_optional( GEN_K ) == boost::none ) {
        config.put(GEN_K, nGens);
    } else {
        nGens = config.get< unsigned int >( GEN_K, nGens);
    }
    std::cerr << "Generations: " << nGens << std::endl;

    population_growth_type pop_grow;
    if( print_config_only ) {
        population_growth_toolkit::getInstance()->tool_configurations( config );

        nGens = 0;
    } else {
        pop_grow = population_growth_toolkit::getInstance()->get_tool( config )->generate();
    }

    simulation_log log( config );

    unsigned int p_size = 0;
    unsigned int i = 0;

    while( i < nGens ) {
        p_size = (*pop_grow)( p_size, i++ );

        sim_engine.simulate(p_size);

        if( log( &sim_engine ) ) {
            log.write();
        }
    }

    boost::property_tree::ptree ofile;
    ofile.add_child( "configuration", config );
    boost::property_tree::write_json( std::cout, ofile );


    return 0;
}
