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
#include "mutation_rate_option.hpp"
#include "config_manager.hpp"

#include <sstream>

const std::string mutation_rate_option::RATE_K = "mu";

mutation_rate_option::mutation_rate_option() {
    config_manager::getInstance()->register_configurable(this);
}

std::string mutation_rate_option::name() const {
    static const std::string n = "Mutation Rate";
    return n;
}

void mutation_rate_option::getOptions( po::options_description & desc ) {
    desc.add_options()
    ( (RATE_K + ",m").c_str(), po::value< mutation_rate_type >()->default_value(0.01), "Mutation Rate")
    ;
}

COMMANDLINE_ACTION mutation_rate_option::validate( const po::variables_map & vm ) {
    return COMMANDLINE_SUCCESS;
}

mutation_rate_option::~mutation_rate_option() {
    config_manager::getInstance()->unregister_configurable( this->name() );
}

static mutation_rate_option c_mut_rate_options;
