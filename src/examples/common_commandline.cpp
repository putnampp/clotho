#include "clotho_commandline.h"

/// GENERAL OPTION KEYS
const string HELP_K = "help";
const string VERSION_K = "version";

/// SIMULATION OPTION KEYS
const string FOUNDER_SIZE_K = "founder-size";
const string MUTATION_RATE_K = "mu";
const string RECOMBINATION_RATE_K = "rho";

/// I/O OPTION KEYS
const string OUTPUT_K = "output";

int parse_commandline( int argc, char ** argv, po::variables_map & vm ) {
    po::options_description gen( "General" );
    gen.add_options()
    ( (HELP_K + ",h").c_str(), "Print this" )
    ( (VERSION_K + ",v").c_str(), "Version" )
    ;

    po::options_description clotho_app( "Simulation Parameters" );
    clotho_app.add_options()
    ( GENERATIONS_K.c_str(), po::value<unsigned int>()->default_value(-1), "Simulate a number of generations.")
    ( FOUNDER_SIZE_K.c_str(), po::value< unsigned int >()->default_value(10000), "Founding population size" )
    ( MUTATION_RATE_K.c_str(), po::value< double >()->default_value( 0.0001), "Mutation rate" )
    ( RECOMBINATION_RATE_K.c_str(), po::value< double>()->default_value( 0.0001 ), "Recombination rate" )
    ;

    po::options_description io_param("I/O parameters");
    io_param.add_options()
    ( (OUTPUT_K + ",o").c_str(), po::value< std::string >()->default_value(""), "Prefix for simulation output files")
    ;

    po::options_description cmdline;

    cmdline.add(gen).add(simulation).add( clotho_app ).add(io_param);
    po::store( po::command_line_parser( argc, argv ).options( cmdline ).run(), vm );

    int res = 0;
    if( vm.count( HELP_K ) ) {
        std::cout << cmdline << std::endl;
        res = 1;
    }

    return res;
}
