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
#include "generations_option.hpp"
#include "config_manager.hpp"

#include <sstream>

const std::string generations_option::SIZE_K = "generations";

generations_option::generations_option() {
    config_manager::getInstance()->register_configurable(this);
}

std::string generations_option::name() const {
    static const std::string n = "Generations";
    return n;
}

void generations_option::getOptions( po::options_description & desc ) {
    desc.add_options()
    ( (SIZE_K + ",g").c_str(), po::value< generations_type >()->default_value(10000), "Number of population generations to simulate")
    ;
}

COMMANDLINE_ACTION generations_option::validate( const po::variables_map & vm ) {
    return COMMANDLINE_SUCCESS;
}

generations_option::~generations_option() {
    config_manager::getInstance()->unregister_configurable( this->name() );
}

static generations_option c_pop_size_option;
