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

#include "clotho/genetics/predefined_population_growth.hpp"
#include <cassert>

const std::string PREDEF_POP_NAME = "predefined";

predefined_population_growth::predefined_population_growth( const std::vector< unsigned int > & sizes ) :
    m_gen_size( sizes )
{
    assert( !m_gen_size.empty() );
}

unsigned int predefined_population_growth::operator()( unsigned int, unsigned int generation ) {
    return m_gen_size[ generation % m_gen_size.size() ];
}

const std::string predefined_population_growth::name() const {
    return PREDEF_POP_NAME;
}

void predefined_population_growth::log( std::ostream & out ) const {
    out << "{" << PREDEF_POP_NAME;
    for(std::vector< unsigned int >::const_iterator cit = m_gen_size.begin(); cit != m_gen_size.end(); ++cit ) {
        out << "," << *cit;
    }
    
    out << "}";
}

predefined_population_growth::~predefined_population_growth() {}
