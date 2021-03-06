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
#include "common_commandline.h"

#include <cassert>
#include <fstream>
#include <boost/algorithm/string.hpp>

#include "simulation_config.h"

#include "clotho/mutation/mutation_rate_parameter.hpp"
#include "clotho/recombination/recombination_rate_parameter.hpp"
#include "clotho/random/seed_parameter.hpp"
#include "logging_parameter.hpp"
#include "population_parameter.hpp"
#include "generation_parameter.hpp"

/// GENERAL OPTION KEYS
const string HELP_K = "help";
const string VERSION_K = "version";

/// SIMULATION OPTION KEYS
//const string GENERATIONS_K = "generations";
const string FOUNDER_SIZE_K = "founder-size";
//const string MUTATION_RATE_K = "mu";
//const string RECOMBINATION_RATE_K = "rho";
//const string REPEAT_K = "repeat";

//const string RNG_SEED_K = "seed";
const string LOG_PERIOD_K = "log-period";

/// I/O OPTION KEYS
const string OUTPUT_K = "output";
const string CONFIG_K = "config";
const string PREFIX_K = "prefix";

int parse_commandline( int argc, char ** argv, po::variables_map & vm ) {
    po::options_description gen( "General" );
    gen.add_options()
    ( (HELP_K + ",h").c_str(), "Print this" )
    ( (VERSION_K + ",v").c_str(), "Version" )
    ;

    typedef double real_type;
    typedef mutation_rate_parameter< real_type > mutation_type;
    typedef recombination_rate_parameter< real_type > recombination_type;

    po::options_description clotho_app( "Simulation Parameters" );
    clotho_app.add_options()
    ( (GENERATIONS_K + ",g").c_str(), po::value<unsigned int>()->default_value(0), "A positive number of generations. A value 0 will return a default configuration file")
    ( (FOUNDER_SIZE_K+ ",P").c_str(), po::value< unsigned int >()->default_value(population_parameter::DEFAULT_POPULATION_SIZE), "Founding population size" )
    ( (MU_K + ",m").c_str(), po::value< real_type >()->default_value( mutation_type::DEFAULT_MUTATION_RATE ), "Mutation rate" )
    ( (RHO_K + ",r").c_str(), po::value< real_type >()->default_value( recombination_type::DEFAULT_RECOMB_RATE ), "Recombination rate" )
    ( (REPETITION_K + ",R").c_str(), po::value< unsigned int >()->default_value( 1 ), "Repetitions" )
    ( (SEED_K + ",s").c_str(), po::value< typename simulation_config::seed_type >()->default_value( seed_parameter<>::DEFAULT_SEED ), "Random number generator initial seed value" )
    ( (PERIOD_K + ",l").c_str(), po::value< unsigned int >()->default_value( -1 ), "Number of generations between population stats calculations. (-1 => only final population)" )
    ;

    po::options_description io_param("I/O parameters");
    io_param.add_options()
    ( (PREFIX_K + ",p").c_str(), po::value< std::string >()->default_value(""), "Prefix for simulation output files")
    ( (CONFIG_K + ",c").c_str(), po::value< std::string >()->default_value(""), "Simulation configuration file")
    ;

    po::options_description cmdline;

    cmdline.add(gen).add( clotho_app ).add(io_param);
    po::store( po::command_line_parser( argc, argv ).options( cmdline ).run(), vm );

    int res = 0;
    if( vm.count( HELP_K ) ) {
        std::cout << cmdline << std::endl;
        res = 1;
    }

    return res;
}

int parse_commandline( int argc, char ** argv, simulation_config & cfg ) {
    po::variables_map vm;

    int res = parse_commandline( argc, argv, vm );
    if( res ) {
        return res;
    }

    cfg.out_path = vm[ PREFIX_K ].as< string >();
    cfg.cfg_path = vm[ CONFIG_K ].as<string>();

    if( !cfg.cfg_path.empty() ) {
        parse_config( cfg );
    } else {
        update_config( vm, cfg );
    }

    return res;
}

int parse_commandline( int argc, char ** argv, boost::property_tree::ptree & cfg ) {
    po::variables_map vm;

    int res = parse_commandline( argc, argv, vm );
    if( res ) {
        return res;
    }

    string cfg_path = vm[ CONFIG_K ].as<string>();
    string out_path = vm[ PREFIX_K ].as<string>();

    if( !cfg_path.empty() ) {
        assert( boost::algorithm::iends_with( cfg_path, ".json") );
        boost::property_tree::read_json(cfg_path, cfg);

        cfg.put( CONFIG_BLOCK_K + "." + CONFIG_K, cfg_path );
        cfg.put( CONFIG_BLOCK_K + "." + PREFIX_K, out_path );
    } else {
        update_config( vm, cfg );
    }

    return res;
}

void update_config( po::variables_map & vm, simulation_config & cfg ) {
    typedef double real_type;

    cfg.nGen = vm[ GENERATIONS_K ].as< unsigned int >();
    cfg.nPop = vm[ FOUNDER_SIZE_K ].as< unsigned int >();
    cfg.nRep = vm[ REPETITION_K ].as< unsigned int >();

    cfg.mu = vm[ MU_K ].as< real_type >();
    cfg.rho = vm[ RHO_K ].as< real_type >();

    cfg.seed = vm[ SEED_K ].as< unsigned int >();
    cfg.log_period = vm[ PERIOD_K ].as< unsigned int >();

    cfg.cfg_path = vm[ CONFIG_K ].as<string>();
    cfg.out_path = vm[ PREFIX_K ].as<string>();
}

void update_config( po::variables_map & vm, boost::property_tree::ptree & cfg ) {
    simulation_config tmp;
    update_config( vm, tmp );

//    if( cfg.get_child_optional( CONFIG_BLOCK_K ) != boost::none ) {
//        add_config(cfg.get_child( CONFIG_BLOCK_K ), tmp);
//    } else {
//        boost::property_tree::ptree t;
//        add_config( t, tmp );
//        cfg.put_child( CONFIG_BLOCK_K, t );
//    }
//
    boost::property_tree::ptree t;
    t = cfg.get_child( CONFIG_BLOCK_K, t );
    add_config( t, tmp );

    cfg.put_child( CONFIG_BLOCK_K, t );
}
