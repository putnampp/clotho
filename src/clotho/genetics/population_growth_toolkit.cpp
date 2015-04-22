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
#include "clotho/genetics/population_growth_toolkit.hpp"

const std::string POPULATION_GROWTH_BLOCK_K = "population_growth";

population_growth_toolkit::population_growth_toolkit() {}

void population_growth_toolkit::tool_configurations( boost::property_tree::ptree & config ) {

    boost::property_tree::ptree t;
    for( generator_citerator it = m_tools.begin(); it != m_tools.end(); ++it ) {
        boost::property_tree::ptree c;
        std::shared_ptr< ipopulation_growth_generator > tmp( it->second->create(c) );

        boost::property_tree::ptree d;
        d.put( "name", it->first );
        d.put_child( "params", c );

        t.put_child( it->first, d );
    }
    config.put_child( POPULATION_GROWTH_BLOCK_K + ".toolkit", t );
}

std::shared_ptr< ipopulation_growth_generator > population_growth_toolkit::get_tool( boost::property_tree::ptree & config ) {
    if( config.get_child_optional( POPULATION_GROWTH_BLOCK_K ) == boost::none ) {
        boost::property_tree::ptree tmp;
        tmp.put( "name", "" );
        tmp.put( "params", "" );
        config.put_child( POPULATION_GROWTH_BLOCK_K, tmp );
        return std::shared_ptr< ipopulation_growth_generator >();
    }

    std::string tname = config.get< std::string >( POPULATION_GROWTH_BLOCK_K + ".name", "" );

    std::shared_ptr< ipopulation_growth_generator > ret;
    if( !tname.empty() ) {
        generator_iterator it = m_tools.find(tname);
        if( it != m_tools.end() ) {
            if( config.get_child_optional( POPULATION_GROWTH_BLOCK_K + ".params" ) == boost::none ) {
                boost::property_tree::ptree t;
                ret = it->second->create( t );
                config.put_child( POPULATION_GROWTH_BLOCK_K + ".params", t );
            } else {
                ret = it->second->create( config.get_child( POPULATION_GROWTH_BLOCK_K + ".params" ) );
            }
        }
    }

    return ret;
}

void population_growth_toolkit::register_tool( ipopulation_growth_generator * gen ) {
    if( gen == NULL ) return;

    generator_iterator it = m_tools.find( gen->name() );
    if( it == m_tools.end() ) {
        m_tools.insert( std::make_pair( gen->name(), gen ) );
    }
}

population_growth_toolkit::~population_growth_toolkit() {}

std::ostream & operator<<( std::ostream & out, const population_growth_toolkit & ftk ) {
    for( population_growth_toolkit::generator_citerator it = ftk.m_tools.begin(); it != ftk.m_tools.end(); ++it ) {
        out << it->first << std::endl;
    }
    return out;
}
