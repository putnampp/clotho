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
#ifndef SIMULATION_LOG_HPP_
#define SIMULATION_LOG_HPP_

#include <boost/property_tree/ptree.hpp>
#include <cstring>
#include <sstream>

#include "state_object.hpp"

class simulation_log {
public:
    typedef boost::property_tree::ptree record_type;

    simulation_log( boost::property_tree::ptree & config );

    void set_path_prefix( std::string & prefix );
    std::string get_path_prefix() const;

    bool operator()( iStateObject * obj );

    void add_record( std::string name, const record_type & rec );

    void purge();

    template < class FillType >
    std::string make_path( FillType suffix ) {
        std::ostringstream oss;

        oss << m_prefix << "_" << suffix << ".json";

        return oss.str();
    }

    void write( std::ostream & out );
    void write();
    void write( std::string & path );

    virtual ~simulation_log();
protected:
    void parse_configuration( boost::property_tree::ptree & config );

    bool m_activated;
    std::string m_prefix;

    boost::property_tree::ptree m_log;
    unsigned int m_count, m_frequency;
    unsigned int m_log_idx;
};

#endif  // SIMULATION_LOG_HPP_
