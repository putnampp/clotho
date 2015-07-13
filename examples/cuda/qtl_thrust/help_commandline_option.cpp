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
#include "help_commandline_option.hpp"

#include "config_manager.hpp"

const std::string HELP_K = "help";

struct help_commandline_option : public iconfigurable {
    help_commandline_option() {
        config_manager::getInstance()->register_configurable( this );
    }

    std::string name() const {
        static std::string n = "Help";
        return n;
    }

    void getOptions( po::options_description & desc ) {
        desc.add_options()
        ( (HELP_K + ",h").c_str(), "Print this" )
        ;
    }

    COMMANDLINE_ACTION validate( const po::variables_map & vm ) {
        if( vm.count( HELP_K ) ) {
            return PRINT_COMMANDLINE;
        }
        return COMMANDLINE_SUCCESS;
    }

    virtual ~help_commandline_option() {
        config_manager::getInstance()->unregister_configurable( this->name() );
    }
};

static help_commandline_option opt_help;
