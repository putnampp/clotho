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
#include "log_writer.hpp"

void OutputLogWriter::write( boost::property_tree::ptree & plog ) {
    boost::property_tree::write_json( std::cout, plog);
}

JSONLogWriter::JSONLogWriter( std::string prefix, std::string suffix ) : 
    m_prefix( prefix )
    , m_suffix( suffix )
    , m_index( 0 )
{}

void JSONLogWriter::write( boost::property_tree::ptree & log ) {
    std::ostringstream oss;

    oss << m_prefix << "." << m_index++ << m_suffix << ".json";
    boost::property_tree::write_json( oss.str(), log);
}
