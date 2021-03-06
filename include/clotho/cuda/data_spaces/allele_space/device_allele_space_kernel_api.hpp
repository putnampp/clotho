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

template < class RealType >
__device__ bool  _resize_space_impl( device_allele_space< RealType > * aspace, unsigned int N );

/**
 *
 *
 *
 */
template < class RealType >
__global__ void _update_free_count( device_allele_space< RealType > * aspace );

/**
 *
 *
 */
template < class RealType >
__global__ void _generate_allele_space( device_allele_space< RealType > * aspace );

/**
 * Merge of allele spaces takes two spaces
 *
 * The output allele space will contain the elements from A and B spaces
 * Elements will be interlaced according to their free state
 *  - If an element is indicated as free in A, then it will be replaced by a non-free element of B.
 *  - Elements will be merged according to the defined order (OrderedTag) of the space
 */
template < class RealType >
__global__ void _merge_allele_space( device_allele_space< RealType > * aspace
                                    , device_allele_space< RealType > * bspace
                                    , device_allele_space< RealType > * output );

#include "clotho/cuda/data_spaces/free_space/device_free_space.hpp"
#include "clotho/cuda/data_spaces/event_space/device_event_space.hpp"

template < class RealType, class IntType, class OrderTag >
__global__ void _merge_space( device_allele_space< RealType > * in_space
                            , device_free_space< IntType, OrderTag > * fspace
                            , device_event_space< IntType, OrderTag > * evts
                            , device_allele_space< RealType > * out_space );

template < class AlleleSpaceType >
__device__ bool _move_allele( AlleleSpaceType * in_space, unsigned int in_idx, AlleleSpaceType * out_space, unsigned int out_idx );

#endif  // DEVICE_ALLELE_SPACE_KERNEL_API_HPP_
