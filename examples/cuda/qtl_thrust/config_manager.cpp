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
#include "config_manager.hpp"

#include <cassert>
#include <iostream>

config_manager::config_manager() {}

config_manager::~config_manager() {}

config_manager * config_manager::getInstance() {
    static config_manager inst;
    return &inst;
}

void config_manager::register_configurable( iconfigurable * obj ) {
    object_iterator it = m_configs.find( obj->name() );
    if( it != m_configs.end() ) {
        std::cerr << "ERROR: Unable to register object configuration.  Configurable object " << obj->name() << " name conflicts\n";
        assert(false);
    }
    m_configs.insert( std::make_pair(obj->name(), obj));
}

void config_manager::unregister_configurable( const std::string & name ) {
    object_iterator it = m_configs.find( name );
    if( it != m_configs.end() ) {
        m_configs.erase( it );
    }
}

int config_manager::parse_commandline( int argc, char ** argv ) {
    po::variables_map vm;
    return parse_commandline( argc, argv, vm );
}

int config_manager::parse_commandline( int argc, char ** argv, po::variables_map & vm ) {
    po::options_description cmdline;
    for( object_iterator it = m_configs.begin(); it != m_configs.end(); ++it ) {
        po::options_description desc( it->first );

        it->second->getOptions( desc );
        cmdline.add( desc );
    }

    po::store( po::command_line_parser( argc, argv ).options( cmdline ).run(), vm );

    int valid = 0;

    for( object_iterator it = m_configs.begin(); it != m_configs.end(); ++it ) {
        COMMANDLINE_ACTION res = it->second->validate( vm );
        if( res == PRINT_MESSAGE ) {
            it->second->printMessage( std::cout );
        } else if( res == PRINT_WARNING ) {
            std::cerr << "WARNING: ";
            it->second->printWarning( std::cerr );
        } else if( res == PRINT_ERROR ) {
            std::cerr << "ERROR: ";
            valid = -2;
            break;
        } else if( res == PRINT_COMMANDLINE ) {
            std::cout << cmdline << "\n";

            valid = -1;
            break;
        }
    }

    return valid;
}
