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
#ifndef CLOTHO_POPULATION_SPACE_ROW_HPP_
#define CLOTHO_POPULATION_SPACE_ROW_HPP_

#include "clotho/data_spaces/phenotype_evaluator/trait_space_vector.hpp"

#include "clotho/utility/bit_helper.hpp"
#include <vector>

namespace clotho {
namespace genetics {


#ifndef USE_VECTOR_ARRAY

#include "clotho/data_spaces/population_space/population_space_row_block.hpp"

#else // USE_VECTOR_ARRAY
#include "clotho/data_spaces/population_space/population_space_row_vector_array.hpp"
#endif  // USE_VECTOR_ARRAY


}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_POPULATION_SPACE_ROW_HPP_

