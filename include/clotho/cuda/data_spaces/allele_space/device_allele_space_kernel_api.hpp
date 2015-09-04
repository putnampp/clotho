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
#ifndef DEVICE_ALLELE_SPACE_KERNEL_API_HPP_
#define DEVICE_ALLELE_SPACE_KERNEL_API_HPP_

#include "clotho/cuda/data_spaces/allele_space/device_allele_space_def.hpp"
#include "clotho/cuda/data_spaces/data_space_kernel_api.hpp"

template < class RealType, class IntType, class OrderTag >
__device__ void _resize_space_impl( device_allele_space< RealType, IntType, OrderTag > * aspace, unsigned int N );

/**
 *
 *
 *
 */
template < class RealType, class IntType, class OrderedTag >
__global__ void _update_free_count( device_allele_space< RealType, IntType, OrderedTag > * aspace );

/**
 *
 *
 */
template < class RealType, class IntType, class OrderedTag >
__global__ void _generate_allele_space( device_allele_space< RealType, IntType, OrderedTag > * aspace );

/**
 * Merge of allele spaces takes two spaces
 *
 * The output allele space will contain the elements from A and B spaces
 * Elements will be interlaced according to their free state
 *  - If an element is indicated as free in A, then it will be replaced by a non-free element of B.
 *  - Elements will be merged according to the defined order (OrderedTag) of the space
 */
template < class RealType, class IntType, class OrderedTag >
__global__ void _merge_allele_space( device_allele_space< RealType, IntType, OrderedTag > * aspace
                                    , device_allele_space< RealType, IntType, OrderedTag > * bspace
                                    , device_allele_space< RealType, IntType, OrderedTag > * output );

#include "clotho/cuda/data_spaces/event_space/device_event_space.hpp"

template < class RealType, class IntType, class OrderTag >
__global__ void _merge_space( device_allele_space< RealType, IntType, OrderTag > * in_space
                            , device_event_space< IntType, OrderTag > * evts
                            , device_allele_space< RealType, IntType, OrderTag > * out_space );


#endif  // DEVICE_ALLELE_SPACE_KERNEL_API_HPP_
