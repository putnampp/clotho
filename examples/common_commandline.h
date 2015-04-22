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
#ifndef CLOTHO_COMMANDLINE_H_
#define CLOTHO_COMMANDLINE_H_

#include "config.h"

#include <boost/program_options.hpp>
#include "simulation_config.h"

namespace po=boost::program_options;

/// GENERAL OPTION KEYS
extern const string HELP_K;
extern const string VERSION_K;

/// SIMULATION OPTION KEYS
extern const string GENERATIONS_K;
extern const string FOUNDER_SIZE_K;
extern const string MUTATION_RATE_K;
extern const string RECOMBINATION_RATE_K;

extern const string REPEAT_K;

extern const string RNG_SEED_K;
extern const string LOG_PERIOD_K;

/// I/O OPTION KEYS
extern const string OUTPUT_K;
extern const string CONFIG_K;

int parse_commandline( int argc, char ** argv, po::variables_map & vm );
int parse_commandline( int argc, char ** argv, simulation_config & cfg );
int parse_commandline( int argc, char ** argv, boost::property_tree::ptree & cfg );

void update_config( po::variables_map & vm, simulation_config & cfg );
void update_config( po::variables_map & vm, boost::property_tree::ptree & cfg );

#endif  // CLOTHO_COMMANDLINE_H_
