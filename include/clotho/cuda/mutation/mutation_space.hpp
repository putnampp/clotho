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
#ifndef MUTATION_SPACE_HPP_
#define MUTATION_SPACE_HPP_

#include "clotho/cuda/allele_space.hpp"

template < class RealType, class IntType >
struct mutation_space {
    RealType *  location;
    IntType *   free_list_offset;

    unsigned int N;
};

#endif  // MUTATION_SPACE_HPP_
