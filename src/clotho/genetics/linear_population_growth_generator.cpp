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
#include "clotho/genetics/linear_population_growth_generator.hpp"
#include "clotho/genetics/population_growth_toolkit.hpp"

linear_population_growth_generator::linear_population_growth_generator() :
    m_A(1.0)
    , m_B(0.0) 
{
    population_growth_toolkit::getInstance()->register_tool(this);
}

linear_population_growth_generator::linear_population_growth_generator( boost::property_tree::ptree & p ) :
    m_A(1.0)
    , m_B(0.0)
{
    parseConfig( p );
}

std::shared_ptr< ipopulation_growth_generator > linear_population_growth_generator::create( boost::property_tree::ptree & c ) const {
    std::shared_ptr< ipopulation_growth_generator > t( new linear_population_growth_generator( c ) );
    return t;
}

std::shared_ptr< ipopulation_growth > linear_population_growth_generator::generate() const {
    std::shared_ptr< ipopulation_growth > t( new linear_population_growth( m_A, m_B ) );
    return t;
}

const std::string linear_population_growth_generator::name() const {
    return LINEAR_POP_NAME;
}

void linear_population_growth_generator::log( std::ostream & out ) const {
    out << "{" << LINEAR_POP_NAME << "_gen"
        << "," << m_A
        << "," << m_B
        << "}";
}

void linear_population_growth_generator::parseConfig( boost::property_tree::ptree & config ) {

    if( config.get_child_optional( "A" ) == boost::none ) {
        config.put( "A", m_A);
    } else {
        m_A = config.get<double>( "A", m_A );
    }

    if( config.get_child_optional( "B" ) == boost::none ) {
        config.put("B", m_B );
    } else {
        m_B = config.get<double>( "B", m_B );
    }
}

linear_population_growth_generator::~linear_population_growth_generator() {}

static const linear_population_growth_generator lpg_gen;
