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
#include "simulation_log.hpp"

#include <iostream>

simulation_log::simulation_log( boost::property_tree::ptree & config ) :
    m_activated( false )
{
    parse_configuration( config );
}

void simulation_log::parse_configuration( boost::property_tree::ptree & config ) {
    boost::property_tree::ptree lconfig;
    if( config.get_child_optional( "logging" ) != boost::none ) {
        lconfig = config.get_child( "logging" );
    }

    if( lconfig.get_child_optional( "activated" ) == boost::none ) {
        lconfig.put( "activated", m_activated );
    } else {
        m_activated = lconfig.get<bool>( "activated", m_activated );
    }

    config.put_child( "logging", lconfig );
}
