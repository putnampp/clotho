#include "clotho/genetics/fitness_toolkit.hpp"

const std::string FITNESS_BLOCK_K = "fitness";
const std::string TOOLKIT_BLOCK_K = "toolkit";
const std::string METRIC_K = "metric";

fitness_toolkit::fitness_toolkit() {}

void fitness_toolkit::tool_configurations( boost::property_tree::ptree & config ) {
    
    for( generator_citerator it = m_tools.begin(); it != m_tools.end(); ++it ) {
        boost::property_tree::ptree c;
        std::shared_ptr< ifitness_generator > tmp( it->second->create(c) );

        config.put_child( FITNESS_BLOCK_K + "." + TOOLKIT_BLOCK_K + "." + it->first, c );
    }
}

std::shared_ptr< ifitness_generator > fitness_toolkit::get_tool( boost::property_tree::ptree & config ) {
    if( config.get_child_optional( FITNESS_BLOCK_K + "." + METRIC_K ) == boost::none ) {
        config.put( FITNESS_BLOCK_K + "." + METRIC_K, "" );
        return std::shared_ptr< ifitness_generator >();
    }

    std::string tname = config.get< std::string >( FITNESS_BLOCK_K + "."  + METRIC_K, "" );

    if( !tname.empty() ) {
        generator_iterator it = m_tools.find(tname);
        if( it != m_tools.end() ) {
            std::string fname = FITNESS_BLOCK_K + "." + TOOLKIT_BLOCK_K + "." + tname;

            boost::property_tree::ptree c;
            if( config.get_child_optional( fname ) != boost::none ) {
                c = config.get_child(fname);
            }

            return it->second->create( c );
        }
    }

    return std::shared_ptr< ifitness_generator >();
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
