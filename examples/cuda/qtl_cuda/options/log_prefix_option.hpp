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
#ifndef LOG_PREFIX_OPTION_HPP_
#define LOG_PREFIX_OPTION_HPP_

#include <string>
#include "clotho/configuration_manager/iconfigurable.hpp"

struct log_prefix_option : public clotho::configuration_manager::iconfigurable {
    typedef std::string path_type;

    static const std::string PREFIX_K;

    log_prefix_option();

    std::string name() const;

    void getOptions( po::options_description & desc );

    clotho::configuration_manager::COMMANDLINE_ACTION validate( const po::variables_map & vm );

    virtual ~log_prefix_option();
};

#endif  // LOG_PREFIX_OPTION_HPP_
