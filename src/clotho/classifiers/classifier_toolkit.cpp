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
#include "clotho/classifiers/classifier_toolkit.hpp"

const std::string CLASSIFIER_BLOCK_K = "classifier";
const std::string TOOLKIT_BLOCK_K = "toolkit";
const std::string TYPE_K = "type";

classifier_toolkit::classifier_toolkit() {}

void classifier_toolkit::tool_configurations( boost::property_tree::ptree & config ) {
    boost::property_tree::ptree t;
    for( generator_citerator it = m_tools.begin(); it != m_tools.end(); ++it ) {
        boost::property_tree::ptree c;

        std::shared_ptr< iclassifier_generator > tmp( it->second->create(c) );
        t.put_child( it->first, c );
    }

    config.put_child( CLASSIFIER_BLOCK_K + "." + TOOLKIT_BLOCK_K, t );
}

std::shared_ptr< iclassifier_generator > classifier_toolkit::get_classifier( boost::property_tree::ptree & config ) {
    if( configs.get_child_optional( CLASSIFIER_BLOCK_K + "." + TYPE_K ) == boost::none ) {
        config.put( CLASSIFIER_BLOCK_K + "." + TYPE_K, "" );
        return std::shared_ptr< iclassifier_generator >();
    }

    std::string tname = config.get< std::string >( CLASSIFIER_BLOCK_K + "." + TYPE_K, "" );

    generator_iterator it = m_tools.find(tname);
    if( it != m_tools.end() ) {
        std::string fname = CLASSIFIER_BLOCK_K + "." + TOOLKIT_BLOCK_K + "." + tname;

        boost::property_tree::ptree c;
        if( config.get_child_optional( fname ) != boost::none ) {
            c = config.get_child(fname);
        }

        return it->second->create( c );
    }

    return std::shared_ptr< iclassifier_generator >();
}

void classifier_toolkit::register_classifier( iclassifier_generator * gen ) {
    if( gen == NULL ) return;

    generator_iterator it = m_classifiers.find( gen->name() );
    if( it == m_classifiers.end() ) {
        m_classifiers.insert( std::make_pair( gen->name(), gen ) );
    }
}

classifier_toolkit::~classifier_toolkit() {}
