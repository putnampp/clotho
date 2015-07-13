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
#include "random_number_options.hpp"
#include "config_manager.hpp"

#include <sstream>

const std::string random_number_options::SEED_K = "seed";

random_number_options::random_number_options() {
    config_manager::getInstance()->register_configurable(this);
}

std::string random_number_options::name() const {
    static const std::string n = "Random Numbers";
    return n;
}

void random_number_options::getOptions( po::options_description & desc ) {
    desc.add_options()
    ( (SEED_K + ",s").c_str(), po::value< seed_type >()->default_value(1234), "Random number generator seed value")
    ;
}

COMMANDLINE_ACTION random_number_options::validate( const po::variables_map & vm ) {
    return COMMANDLINE_SUCCESS;
}

random_number_options::~random_number_options() {
    config_manager::getInstance()->unregister_configurable( this->name() );
}

static random_number_options c_rand_options;
