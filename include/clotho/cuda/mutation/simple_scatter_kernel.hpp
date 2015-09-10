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
#ifndef SCATTER_KERNELS_HPP_
#define SCATTER_KERNELS_HPP_

#include <cuda.h>

template < class StateType, class AlleleSpaceType, class EventSpaceType, class SequenceSpaceType >
__global__ void _scatter_mutation( StateType * states
                                    , AlleleSpaceType   * alleles
                                    , EventSpaceType    * events
                                    , SequenceSpaceType * sequences );

#include "clotho/cuda/mutation/scatter_unordered_impl.hpp"

#endif  // SCATTER_KERNELS_HPP_
