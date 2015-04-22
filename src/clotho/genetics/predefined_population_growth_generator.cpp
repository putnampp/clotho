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
#include "clotho/genetics/predefined_population_growth_generator.hpp"
#include "clotho/genetics/population_growth_toolkit.hpp"

#include "clotho/utility/log_helper.hpp"
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

const unsigned int DEF_POP_SIZE = 10000;

predefined_population_growth_generator::predefined_population_growth_generator() :
    m_gens()
{
    m_gens.push_back(DEF_POP_SIZE);
    population_growth_toolkit::getInstance()->register_tool(this);
}

predefined_population_growth_generator::predefined_population_growth_generator( boost::property_tree::ptree & p ) :
    m_gens()
{
    parseConfig( p );
}

std::shared_ptr< ipopulation_growth_generator > predefined_population_growth_generator::create( boost::property_tree::ptree & c ) const {
    std::shared_ptr< ipopulation_growth_generator > t( new predefined_population_growth_generator( c ) );
    return t;
}

std::shared_ptr< ipopulation_growth > predefined_population_growth_generator::generate() const {
    std::shared_ptr< ipopulation_growth > t( new predefined_population_growth( m_gens ) );
    return t;
}

const std::string predefined_population_growth_generator::name() const {
    return PREDEF_POP_NAME;
}

void predefined_population_growth_generator::log( std::ostream & out ) const {
    out << "{" << PREDEF_POP_NAME << "_gen";

    for( std::vector< unsigned int >::const_iterator cit = m_gens.begin(); cit != m_gens.end(); ++cit ) {
        out << "," << *cit;
    }

    out << "}";
}

void predefined_population_growth_generator::parseConfig( boost::property_tree::ptree & config ) {
    if( config.get_child_optional( "sizes" ) == boost::none ) {
        if( m_gens.empty() ) {
            m_gens.push_back(DEF_POP_SIZE);
        }
        boost::property_tree::ptree t;
        BOOST_FOREACH( auto& v, m_gens ) {
            clotho::utility::add_value_array( t, v );
        }
        config.add_child( "sizes", t );
    } else {
        BOOST_FOREACH( auto & v, config.get_child( "sizes" ) ) {
            if( v.second.empty() ) {    // node should be empty (no children, just data)
                unsigned int s = v.second.get< unsigned int >("", 0 );
                m_gens.push_back(s);
            }
        }
        if( m_gens.empty() ) {
            m_gens.push_back( DEF_POP_SIZE );
        }
    }
}

predefined_population_growth_generator::~predefined_population_growth_generator() {}

static const predefined_population_growth_generator ppg_gen;
