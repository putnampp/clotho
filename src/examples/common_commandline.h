#ifndef CLOTHO_COMMANDLINE_H_
#define CLOTHO_COMMANDLINE_H_

#include "clotho.h"
#include <boost/program_options.hpp>

namespace po=boost::program_options;

/// GENERAL OPTION KEYS
extern const string HELP_K;
extern const string VERSION_K;

/// SIMULATION OPTION KEYS
extern const string GENERATIONS_K;
extern const string FOUNDER_SIZE_K;
extern const string MUTATION_RATE_K;
extern const string RECOMBINATION_RATE_K;

/// I/O OPTION KEYS
extern const string OUTPUT_K;

int parse_commandline( int argc, char ** argv, po::variables_map & vm );

#endif  // CLOTHO_COMMANDLINE_H_
