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
#ifndef CROSSOVER_VALIDATER_HPP_
#define CROSSOVER_VALIDATER_HPP_

#include <boost/property_tree/ptree.hpp>

template < class CrossType >
bool validate( CrossType & ct, boost::property_tree::ptree & err );

#include "validate_crossover_matrix_5.hpp"
#include "validate_crossover_matrix_4.hpp"

#endif  // CROSSOVER_VALIDATER_HPP_
