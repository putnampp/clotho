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
#ifndef GENERATION_PARAMETER_HPP_
#define GENERATION_PARAMETER_HPP_

#include <boost/property_tree/ptree.hpp>
#include <string>
#include "common_strings.h"

extern const std::string GENERATIONS_K;

struct generation_parameter {
    static const unsigned int DEFAULT_GENERATIONS = 10000;

    generation_parameter( unsigned int s = DEFAULT_GENERATIONS );
    generation_parameter( boost::property_tree::ptree & config );

    void write_parameter( boost::property_tree::ptree & l );

    virtual ~generation_parameter();

    unsigned int m_size;
};
#endif  // GENERATION_PARAMETER_HPP_
