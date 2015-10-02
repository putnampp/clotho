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
#ifndef SELECT_AND_CROSSOVER_KERNELS_HPP_
#define SELECT_AND_CROSSOVER_KERNELS_HPP_

#include <cuda.h>

template < class StateType
           , class AlleleSpaceType
           , class SequenceSpaceType
           , class EventSpaceType1
           , class EventSpaceType2
           , class FreeSpaceType >
__global__ void select_and_crossover_kernel( StateType * states
                                            , SequenceSpaceType * parent_seqs
                                            , EventSpaceType1 * parent_ids
                                            , AlleleSpaceType * alleles
                                            , FreeSpaceType * free_space
                                            , EventSpaceType2 * xover_events
                                            , SequenceSpaceType * offspring_seqs );

#include "clotho/cuda/crossover/select_and_crossover_unordered_impl.hpp"

#endif  // SELECT_AND_CROSSOVER_KERNELS_HPP_
