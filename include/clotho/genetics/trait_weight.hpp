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
#ifndef TRAIT_WEIGHT_HPP_
#define TRAIT_WEIGHT_HPP_

#include "clotho/utility/clotho_strings.hpp"
#include <vector>

template < class ValueType = double >
struct trait_weight {
    typedef ValueType                   value_type;
    typedef std::vector< value_type >   vector_type;
};

#include "trait_weight_generator.hpp"

#endif  // TRAIT_WEIGHT_HPP_
