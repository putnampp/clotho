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
#include "population_size_option.hpp"
#include "config_manager.hpp"

#include <sstream>

const std::string population_size_option::SIZE_K = "pop-size";

population_size_option::population_size_option() {
    config_manager::getInstance()->register_configurable(this);
}

std::string population_size_option::name() const {
    static const std::string n = "Population Size";
    return n;
}

void population_size_option::getOptions( po::options_description & desc ) {
    desc.add_options()
    ( (SIZE_K + ",p").c_str(), po::value< population_size_type >()->default_value(0.01), "Population Size")
    ;
}

COMMANDLINE_ACTION population_size_option::validate( const po::variables_map & vm ) {
    return COMMANDLINE_SUCCESS;
}

population_size_option::~population_size_option() {
    config_manager::getInstance()->unregister_configurable( this->name() );
}

static population_size_option c_pop_size_option;
