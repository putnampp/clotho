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
#include "clotho/fitness/fitness_toolkit.hpp"
#include "clotho/utility/clotho_strings.hpp"

#include <boost/foreach.hpp>

fitness_toolkit::fitness_toolkit() {}

void fitness_toolkit::tool_configurations( boost::property_tree::ptree & config ) {

    boost::property_tree::ptree t;
    for( generator_citerator it = m_tools.begin(); it != m_tools.end(); ++it ) {
        boost::property_tree::ptree c;
        std::shared_ptr< ifitness_generator > tmp( it->second->create(c) );

        boost::property_tree::ptree d, p, dep;
        p = c.get_child( FITNESS_BLOCK_K + "." + PARAM_K, p );

        BOOST_FOREACH( auto& v, c ) {
            if( v.first != FITNESS_BLOCK_K ) {
                dep.put_child( v.first, v.second );
            }
        }

        d.put( NAME_K, it->first );

        if( !dep.empty() )
            d.put_child( DEPENDS_K, dep );

        if( !p.empty() )
            d.put_child( PARAM_K, p );

        t.put_child( it->first, d );
    }
    config.put_child( FITNESS_BLOCK_K + ".toolkit", t );
}

std::shared_ptr< ifitness_generator > fitness_toolkit::get_tool( boost::property_tree::ptree & config ) {
    boost::property_tree::ptree fblock;
    fblock = config.get_child( FITNESS_BLOCK_K, fblock );

    std::string tname = fblock.get< std::string >(NAME_K, "");
    std::shared_ptr< ifitness_generator > ret;

    if( fblock.empty() ) {
        fblock.put( NAME_K, "" );
        fblock.put( PARAM_K, "" );
    } else if( !tname.empty() ) {
        generator_iterator it = m_tools.find(tname);
        if( it != m_tools.end() ) {
            ret = it->second->create( config );
        }
    }

    config.put_child( FITNESS_BLOCK_K, fblock );
    return ret;
}

void fitness_toolkit::register_tool( ifitness_generator * gen ) {
    if( gen == NULL ) return;

    generator_iterator it = m_tools.find( gen->name() );
    if( it == m_tools.end() ) {
        m_tools.insert( std::make_pair( gen->name(), gen ) );
    }
}

fitness_toolkit::~fitness_toolkit() {}

std::ostream & operator<<( std::ostream & out, const fitness_toolkit & ftk ) {
    for( fitness_toolkit::generator_citerator it = ftk.m_tools.begin(); it != ftk.m_tools.end(); ++it ) {
        out << it->first << std::endl;
    }
    return out;
}
