#include "common_commandline.h"

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
        cfg.nGen = vm[ GENERATIONS_K ].as< unsigned int >();
        cfg.nPop = vm[ FOUNDER_SIZE_K ].as< unsigned int >();
        cfg.nRep = vm[ REPEAT_K ].as< unsigned int >();

        cfg.mu = vm[ MUTATION_RATE_K ].as< double >();
        cfg.rho = vm[ RECOMBINATION_RATE_K ].as< double >();

        cfg.seed = vm[ RNG_SEED_K ].as< unsigned int >();
        cfg.log_period = vm[ LOG_PERIOD_K ].as< unsigned int >();
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

        cfg.put( CONFIG_BLOCK_K + "." + CONFIG_K, cfg_path);
        cfg.put( CONFIG_BLOCK_K + "." + OUTPUT_K, out_path);
    } else {
        simulation_config tmp_cfg;
        tmp_cfg.cfg_path = cfg_path;
        tmp_cfg.out_path = out_path;
        tmp_cfg.nGen = vm[ GENERATIONS_K ].as< unsigned int >();
        tmp_cfg.nPop = vm[ FOUNDER_SIZE_K ].as< unsigned int >();
        tmp_cfg.nRep = vm[ REPEAT_K ].as< unsigned int >();

        tmp_cfg.mu = vm[ MUTATION_RATE_K ].as< double >();
        tmp_cfg.rho = vm[ RECOMBINATION_RATE_K ].as< double >();

        tmp_cfg.seed = vm[ RNG_SEED_K ].as< unsigned int >();
        tmp_cfg.log_period = vm[ LOG_PERIOD_K ].as< unsigned int >();

        boost::property_tree::ptree t;
        add_config( t, tmp_cfg );
        cfg.add_child( CONFIG_BLOCK_K, t );
    }

    return res;
}
