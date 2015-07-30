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
#ifndef RECOMBINATION_RATE_OPTION_HPP_
#define RECOMBINATION_RATE_OPTION_HPP_

#include <string>

#include "clotho/configuration_manager/iconfigurable.hpp"

struct recombination_rate_option : public clotho::configuration_manager::iconfigurable {
    typedef double recombination_rate_type;

    static const std::string RATE_K;
    recombination_rate_option();
 
    std::string name() const;

    void getOptions( po::options_description & desc );

    clotho::configuration_manager::COMMANDLINE_ACTION validate( const po::variables_map & vm );

    virtual ~recombination_rate_option();
};

#endif  // RECOMBINATION_RATE_OPTION_HPP_
