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
#include "generation_parameter.hpp"

const std::string GENERATIONS_K = "generations";

const unsigned int generation_parameter::DEFAULT_GENERATIONS;

generation_parameter::generation_parameter( unsigned int s ) :
    m_size( s )
{}

generation_parameter::generation_parameter( boost::property_tree::ptree & config ) :
    m_size( DEFAULT_GENERATIONS )
{
    boost::property_tree::ptree lconfig;
    lconfig = config.get_child( GENERATIONS_K, lconfig );

    m_size = lconfig.get< unsigned int >( SIZE_K, m_size );

    lconfig.put( SIZE_K, m_size );
    config.put_child( GENERATIONS_K, lconfig );
}

void generation_parameter::write_parameter( boost::property_tree::ptree & l ) {
    boost::property_tree::ptree c;
    c.put( SIZE_K, m_size );

    l.put_child( GENERATIONS_K, c );
}

generation_parameter::~generation_parameter() {}
