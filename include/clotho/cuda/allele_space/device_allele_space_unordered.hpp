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
#ifndef DEVICE_ALLELE_SPACE_UNORDERED_HPP_
#define DEVICE_ALLELE_SPACE_UNORDERED_HPP_

#include "clotho/cuda/allele_space/device_allele_space_def.hpp"

#include "clotho/cuda/allele_space/tags/unordered_tag.hpp"

template < class RealType, class IntType >
__global__ void _update_free_count( device_allele_space< RealType, IntType, unordered_tag > * aspace ) {

}

#endif  // DEVICE_ALLELE_SPACE_UNORDERED_HPP_
