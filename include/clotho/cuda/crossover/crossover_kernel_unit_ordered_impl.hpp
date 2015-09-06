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
#ifndef CROSSOVER_KERNEL_UNIT_ORDERED_IMPL_HPP_
#define CROSSOVER_KERNEL_UNIT_ORDERED_IMPL_HPP_

#include "clotho/cuda/data_spaces/allele_space/device_allele_space.hpp"
#include "clotho/cuda/data_spaces/sequence_space/device_sequence_space.hpp"
#include "clotho/cuda/distributions/poisson_distribution.hpp"

#include "clotho/cuda/data_spaces/tags/unit_ordered_tag.hpp"

template < class StateType, class RealType, class IntType >
__global__ void crossover_kernel( StateType * states
                                , device_allele_space< RealType, IntType, unit_ordered_tag< IntType > > * alleles
                                , poisson_cdf< RealType, 32 > * events
                                , device_sequence_space< IntType > * sequences ) {

}

#endif  // CROSSOVER_KERNEL_UNIT_ORDERED_IMPL_HPP_
