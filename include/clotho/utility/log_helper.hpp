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
#ifndef LOG_HELPER_HPP_
#define LOG_HELPER_HPP_

#include <boost/property_tree/ptree.hpp>
#include "clotho/utility/timer.hpp"

#include <string>
#include <vector>
#include <memory>

namespace clotho {
namespace utility {

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

template < class Iter >
void add_value_array( boost::property_tree::ptree & array, Iter first, Iter last ) {
    while( first != last ) {
        add_value_array( array, (*first) );
        ++first;
    }
}

template < class A >
void add_value_array( boost::property_tree::ptree & array, const std::vector< A > & t ) {
    boost::property_tree::ptree n;
    add_value_array( n, t.begin(), t.end() );
    array.push_back( std::make_pair("", n ) );
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


#if __cplusplus > 201103L
// shared_ptr only defined in c++11 and higher
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
#endif

}   // namespace utility
}   // namespace clotho

#endif  // LOG_HELPER_HPP_
