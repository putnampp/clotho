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
#include "configuration_option.hpp"
#include "clotho/configuration_manager/config_manager.hpp"

const std::string configuration_option::CONFIG_K = "config";

configuration_option::configuration_option() {
    clotho::configuration_manager::config_manager::getInstance()->register_configurable(this);
}

std::string configuration_option::name() const {
    return CONFIG_K;
}

void configuration_option::getOptions( po::options_description & desc ) {
    desc.add_options()
    ((CONFIG_K + ",c").c_str(), po::value< path_type >(), "Path to JSON formatted configuration file" )
    ;
}

clotho::configuration_manager::COMMANDLINE_ACTION configuration_option::validate( const po::variables_map & vm ) {
    clotho::configuration_manager::COMMANDLINE_ACTION res = clotho::configuration_manager::COMMANDLINE_SUCCESS;

//    if( vm.count( CONFIG_K ) ) {
//        
//    }

    return res;
}

configuration_option::~configuration_option() {
    clotho::configuration_manager::config_manager::getInstance()->unregister_configurable( this->name() );
}

static configuration_option c_conf_option;
