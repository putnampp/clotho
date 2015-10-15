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
#ifndef SEED_PARAMETER_HPP_
#define SEED_PARAMETER_HPP_

#include <boost/property_tree/ptree.hpp>
#include <string>

extern const std::string   RNG_BLOCK_K;
extern const std::string   SEED_K;

struct seed_parameter {
    typedef unsigned long long seed_type;

    static const seed_type DEFAULT_SEED = 0;

    seed_type m_seed;

    seed_parameter( seed_type s = DEFAULT_SEED );
    seed_parameter( boost::property_tree::ptree & config );
 
    void write_parameter( boost::property_tree::ptree & l );

    virtual ~seed_parameter();
};

#endif  // SEED_PARAMETER_HPP_
