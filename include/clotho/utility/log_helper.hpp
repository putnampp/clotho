#ifndef LOG_HELPER_HPP_
#define LOG_HELPER_HPP_

#include <boost/property_tree/ptree.hpp>
#include "clotho/utility/timer.hpp"

#include <string>

template < class Value >
void add_value_array( boost::property_tree::ptree & array, const Value & t ) {
    boost::property_tree::ptree node;

    node.put("", t );
    array.push_back( std::make_pair("", node));
}

template < class A, class B >
void add_value_array( boost::property_tree::ptree & array, const std::pair< A, B > & t ) {
    boost::property_tree::ptree a,b,c;
    a.put( "", t.first );
    b.put( "", t.second );

    c.push_back (std::make_pair( "", a ) );
    c.push_back (std::make_pair( "", b ) );
    array.push_back( std::make_pair("", c ) );
}

inline void add_value_array( boost::property_tree::ptree & array, const clotho::utility::timer & t ) {
    boost::property_tree::ptree node;

    node.put("", t.elapsed().count() );
    array.push_back( std::make_pair("", node));
}

inline void add_node( boost::property_tree::ptree & root, const std::string & path, const boost::property_tree::ptree & node ) {
    if( !node.empty() ) {
        root.add_child( path, node );
    }
}

template < class Sequence >
void add_node( boost::property_tree::ptree & r, const std::string & path, const std::pair< std::shared_ptr< Sequence >, std::shared_ptr< Sequence > > & ind ) {
    std::ostringstream oss;

    oss << ind.first << ": " << (*ind.first);
    r.put( path + ".first", oss.str());

    oss.str("");
    oss.clear();

    oss << ind.second << ": " << (*ind.second);
    r.put( path + ".second", oss.str() );
}

#endif  // LOG_HELPER_HPP_
