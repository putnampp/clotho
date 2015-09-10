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
#include "log_prefix_option.hpp"
#include "clotho/configuration_manager/config_manager.hpp"

const std::string log_prefix_option::PREFIX_K = "prefix";

log_prefix_option::log_prefix_option() {
    clotho::configuration_manager::config_manager::getInstance()->register_configurable(this);
}

std::string log_prefix_option::name() const {
    return PREFIX_K;
}

void log_prefix_option::getOptions( po::options_description & desc ) {
    desc.add_options()
    ((PREFIX_K + ",p").c_str(), po::value< path_type >(), "Prefix for JSON formatted log file" )
    ;
}

clotho::configuration_manager::COMMANDLINE_ACTION log_prefix_option::validate( const po::variables_map & vm ) {
    clotho::configuration_manager::COMMANDLINE_ACTION res = clotho::configuration_manager::COMMANDLINE_SUCCESS;

    return res;
}

log_prefix_option::~log_prefix_option() {
    clotho::configuration_manager::config_manager::getInstance()->unregister_configurable( this->name() );
}

static log_prefix_option c_log_pref_option;
