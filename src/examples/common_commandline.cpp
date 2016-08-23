#include "common_commandline.h"

#include <iostream>
#include <cassert>
#include <fstream>
#include <boost/algorithm/string.hpp>

#include "simulation_config.h"

/// GENERAL OPTION KEYS
const string HELP_K = "help";
const string VERSION_K = "version";

/// SIMULATION OPTION KEYS
const string GENERATIONS_K = "generations";
const string FOUNDER_SIZE_K = "founder-size";
const string MUTATION_RATE_K = "mu";
const string RECOMBINATION_RATE_K = "rho";
const string REPEAT_K = "repeat";
const string RNG_SEED_K = "seed";
const string LOG_PERIOD_K = "log-period";

/// I/O OPTION KEYS
const string OUTPUT_K = "output";
const string CONFIG_K = "config";

int parse_commandline( int argc, char ** argv, po::variables_map & vm ) {
    po::options_description gen( "General" );
    gen.add_options()
    ( (HELP_K + ",h").c_str(), "Print this" )
    ( (VERSION_K + ",v").c_str(), "Version" )
    ;

    po::options_description clotho_app( "Simulation Parameters" );
    clotho_app.add_options()
    ( (GENERATIONS_K + ",g").c_str(), po::value<unsigned int>()->default_value(DEFAULT_GENERATIONS), "Simulate a number of generations.")
    ( (FOUNDER_SIZE_K+ ",p").c_str(), po::value< unsigned int >()->default_value(DEFAULT_POPULATION_SIZE), "Founding population size" )
    ( (MUTATION_RATE_K + ",m").c_str(), po::value< double >()->default_value( DEFAULT_MUTATION_RATE ), "Mutation rate" )
    ( (RECOMBINATION_RATE_K + ",r").c_str(), po::value< double>()->default_value( DEFAULT_RECOMB_RATE ), "Recombination rate" )
    ( (REPEAT_K + ",R").c_str(), po::value< unsigned int >()->default_value( 1 ), "Repetitions" )
    ( (RNG_SEED_K + ",s").c_str(), po::value< unsigned int >()->default_value( DEFAULT_SEED ), "Random number generator initial seed value" )
    ( (LOG_PERIOD_K + ",l").c_str(), po::value< unsigned int >()->default_value( -1 ), "Number of generations between population stats calculations. (-1 => only final population)" )
    ;

    po::options_description io_param("I/O parameters");
    io_param.add_options()
    ( (OUTPUT_K + ",o").c_str(), po::value< std::string >()->default_value(""), "Prefix for simulation output files")
    ( (CONFIG_K + ",i").c_str(), po::value< std::string >()->default_value(""), "Simulation configuration file")
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

    cfg.out_path = vm[ OUTPUT_K ].as< string >();
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
    string out_path = vm[ OUTPUT_K ].as<string>();

    if( !cfg_path.empty() ) {
        assert( boost::algorithm::iends_with( cfg_path, ".json") );
        boost::property_tree::read_json(cfg_path, cfg);

        cfg.put( CONFIG_BLOCK_K + "." + CONFIG_K, cfg_path );
        cfg.put( CONFIG_BLOCK_K + "." + OUTPUT_K, out_path );
    } else {
        update_config( vm, cfg );
    }

    return res;
}

void update_config( po::variables_map & vm, simulation_config & cfg ) {
    cfg.nGen = vm[ GENERATIONS_K ].as< unsigned int >();
    cfg.nPop = vm[ FOUNDER_SIZE_K ].as< unsigned int >();
    cfg.nRep = vm[ REPEAT_K ].as< unsigned int >();

    cfg.mu = vm[ MUTATION_RATE_K ].as< double >();
    cfg.rho = vm[ RECOMBINATION_RATE_K ].as< double >();

    cfg.seed = vm[ RNG_SEED_K ].as< unsigned int >();
    cfg.log_period = vm[ LOG_PERIOD_K ].as< unsigned int >();

    cfg.cfg_path = vm[ CONFIG_K ].as<string>();
    cfg.out_path = vm[ OUTPUT_K ].as<string>();
}

void update_config( po::variables_map & vm, boost::property_tree::ptree & cfg ) {
    simulation_config tmp;
    update_config( vm, tmp );

    if( cfg.get_child_optional( CONFIG_BLOCK_K ) != boost::none ) {
        add_config(cfg.get_child( CONFIG_BLOCK_K ), tmp);
    } else {
        boost::property_tree::ptree t;
        add_config( t, tmp );
        cfg.put_child( CONFIG_BLOCK_K, t );
    }
}
