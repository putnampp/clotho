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
#ifndef RANDOM_NUMBER_OPTIONS_HPP_
#define RANDOM_NUMBER_OPTIONS_HPP_

#include <string>

#include "iconfigurable.hpp"

struct random_number_options : public iconfigurable {
    typedef unsigned long long seed_type;

    static const std::string SEED_K;
    random_number_options();
 
    std::string name() const;

    void getOptions( po::options_description & desc );

    COMMANDLINE_ACTION validate( const po::variables_map & vm );

    virtual ~random_number_options();
};

#endif  // RANDOM_NUMBER_OPTIONS_HPP_
