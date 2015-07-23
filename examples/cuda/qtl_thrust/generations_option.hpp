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
#ifndef GENERATIONS_OPTION_HPP_
#define GENERATIONS_OPTION_HPP_

#include <string>

#include "iconfigurable.hpp"

struct generations_option : public iconfigurable {
    typedef unsigned int generations_type;

    static const std::string SIZE_K;
    generations_option();
 
    std::string name() const;

    void getOptions( po::options_description & desc );

    COMMANDLINE_ACTION validate( const po::variables_map & vm );

    virtual ~generations_option();
};

#endif  // GENERATIONS_OPTION_HPP_
