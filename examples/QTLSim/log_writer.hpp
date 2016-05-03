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
#ifndef LOG_WRITER_HPP_
#define LOG_WRITER_HPP_

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <iostream>
#include <sstream>
#include <memory>
#include <cstring>

struct log_writer {
    virtual void write( boost::property_tree::ptree & ) = 0;
};

struct OutputLogWriter : public log_writer {
    void write( boost::property_tree::ptree & plog ); 
};

struct JSONLogWriter : public log_writer {
    std::string m_prefix, m_suffix;
    unsigned int m_index;

    JSONLogWriter( std::string prefix, std::string suffix );

    void write( boost::property_tree::ptree & log );
};

static std::shared_ptr< log_writer > makeLogWriter( std::string prefix, std::string suffix ) {
    if( prefix.empty() ) {
        std::shared_ptr< log_writer > logger( new OutputLogWriter() );
        return logger;
    } else {
        std::shared_ptr< log_writer > logger( new JSONLogWriter( prefix, suffix ) );
        return logger;
    }
}

#endif  // LOG_WRITER_HPP_
