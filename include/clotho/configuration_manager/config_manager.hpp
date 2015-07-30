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
#ifndef CONFIG_MANAGER_HPP_
#define CONFIG_MANAGER_HPP_

#include <boost/program_options.hpp>

#include "commandline_actions.hpp"
#include "iconfigurable.hpp"

namespace clotho {
namespace configuration_manager {

class config_manager {
public:
    typedef std::map< std::string, iconfigurable * > object_map;
    typedef typename object_map::iterator object_iterator;

    static config_manager * getInstance();
    config_manager();

    void register_configurable( iconfigurable * obj );
    void unregister_configurable( const std::string & n );

    int parse_commandline( int argc, char ** argv );
    int parse_commandline( int argc, char ** argv, po::variables_map & vm );

    virtual ~config_manager();
protected:
    object_map m_configs;
};

}   // namespace configuration_manager
}   // namespace clotho

#endif  // CONFIG_MANAGER_HPP_
