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
#include "recombination_rate_option.hpp"
#include "clotho/configuration_manager/config_manager.hpp"

#include <sstream>

const std::string recombination_rate_option::RATE_K = "rho";

recombination_rate_option::recombination_rate_option() {
    clotho::configuration_manager::config_manager::getInstance()->register_configurable(this);
}

std::string recombination_rate_option::name() const {
    static const std::string n = "Recombination Rate";
    return n;
}

void recombination_rate_option::getOptions( po::options_description & desc ) {
    desc.add_options()
    ( (RATE_K + ",r").c_str(), po::value< recombination_rate_type >()->default_value(0.01), "Recombination Rate")
    ;
}

clotho::configuration_manager::COMMANDLINE_ACTION recombination_rate_option::validate( const po::variables_map & vm ) {
    return clotho::configuration_manager::COMMANDLINE_SUCCESS;
}

recombination_rate_option::~recombination_rate_option() {
    clotho::configuration_manager::config_manager::getInstance()->unregister_configurable( this->name() );
}

static recombination_rate_option c_rec_rate_options;
