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
#include "qtl_thrust_general_options.hpp"

#include "clotho/configuration_manager/config_manager.hpp"

#include <sstream>

const std::string VERSION_K = "version";
const std::string DESCRIPTION_K = "description";

struct qtl_thrust_general_options : public clotho::configuration_manager::iconfigurable {

    qtl_thrust_general_options() {
        clotho::configuration_manager::config_manager::getInstance()->register_configurable( this );
    }

    std::string name() const {
        static const std::string n = "General";
        return n;
    }

    void getOptions( po::options_description & desc ) {
        desc.add_options()
        ( VERSION_K.c_str(), "software version" )
        ( (DESCRIPTION_K).c_str(), "Program Description" )
        ;
    }

    clotho::configuration_manager::COMMANDLINE_ACTION validate( const po::variables_map & vm ) {
        clotho::configuration_manager::COMMANDLINE_ACTION res = clotho::configuration_manager::COMMANDLINE_SUCCESS;

        if( vm.count( VERSION_K ) ) {
            res = clotho::configuration_manager::PRINT_MESSAGE;
            msg << MAJOR_VERSION << "." << MINOR_VERSION << "." << REVISION << "\n"; 
        }

        if( vm.count( DESCRIPTION_K ) ) {
            res = clotho::configuration_manager::PRINT_MESSAGE;

            msg << "This program does something cool\n";
        }

        return res;
    }

    void printMessage( std::ostream & out ) {
        out << msg.str();
        
        // reset message
        msg.str("");
        msg.clear();
    }

    virtual ~qtl_thrust_general_options() {
        clotho::configuration_manager::config_manager::getInstance()->unregister_configurable( this->name());
    }

    std::ostringstream msg;
};

static qtl_thrust_general_options c_qtl_thrust;
