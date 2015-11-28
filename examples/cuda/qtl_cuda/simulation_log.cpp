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

#include <boost/property_tree/json_parser.hpp>
#include <sstream>

simulation_log::simulation_log( boost::property_tree::ptree & config, std::string prefix ) :
    logging_parameter( config )
    , m_activated( false )
    , m_prefix(prefix)
    , m_count(0)
//    , m_frequency(0)
    , m_log_idx(0)
{
//    parse_configuration( config );
    m_count = m_period - 1;
}

void simulation_log::set_path_prefix( std::string & prefix ) {
    m_prefix = prefix;
}

std::string simulation_log::get_path_prefix() const {
    return m_prefix;
}

bool simulation_log::operator()( clotho::utility::iStateObject * obj ) {
    if( m_count-- == 0 ) {
        obj->get_state( m_log );
        m_count = m_period - 1;

        return true;
    }
    return false;
}

void simulation_log::add_record( std::string name, const record_type & rec ) {

    m_log.add_child( name, rec );
}

void simulation_log::write( std::ostream & out ) {
    boost::property_tree::write_json( out, m_log );

    purge();
}

void simulation_log::write() {
    std::string p = make_path( (++m_log_idx) );
    write( p );
}

void simulation_log::write( std::string & path ) {
    boost::property_tree::write_json( path, m_log );
    purge();
}

void simulation_log::purge() {
    m_log.clear();
}

/*
void simulation_log::parse_configuration( boost::property_tree::ptree & config ) {
    boost::property_tree::ptree lconfig;
    if( config.get_child_optional( "logging" ) != boost::none ) {
        lconfig = config.get_child( "logging" );
    }

    if( lconfig.get_child_optional( "path_prefix" ) == boost::none ) {
        lconfig.put( "path_prefix", m_prefix );
    } else {
        m_prefix = lconfig.get< std::string >( "path_prefix", m_prefix );
    }

    if( lconfig.get_child_optional( "frequency" ) == boost::none ) {
        lconfig.put( "frequency", m_frequency );
    } else {
        m_frequency = lconfig.get< unsigned int >( "frequency", m_frequency );

        m_count = m_frequency - 1;
    }

    config.put_child( "logging", lconfig );
}*/

simulation_log::~simulation_log() { }
