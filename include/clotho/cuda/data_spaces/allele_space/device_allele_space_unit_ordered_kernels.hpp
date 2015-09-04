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
#ifndef DEVICE_ALLELE_SPACE_UNIT_ORDERED_KERNELS_HPP_
#define DEVICE_ALLELE_SPACE_UNIT_ORDERED_KERNELS_HPP_

#include "clotho/cuda/data_spaces/allele_space/device_allele_space_unit_ordered.hpp"
#include "clotho/cuda/data_spaces/allele_space/device_allele_space_kernel_api.hpp"

template < class RealType, class IntType >
__global__ void _update_free_count( device_allele_space< RealType, IntType, unit_ordered_tag< IntType > > * aspace ) {
    IntType tid = threadIdx.y *blockDim.x + threadIdx.x;
    IntType lane_mask = (1 << (tid & 31));

    IntType * local_free_list = aspace->free_list;
    unsigned int * fcounts = aspace->free_count;

    unsigned int N = aspace->free_list_size();
    unsigned int _count = 0;

    while( N-- ) {
        IntType f = *local_free_list;

        _count += !!(f & lane_mask);
        ++local_free_list;
    }

    fcounts[tid] = _count;
}


template < class RealType, class IntType >
__global__ void _merge_space( device_allele_space< RealType, IntType, unit_ordered_tag< IntType > > * in_space
                            , device_event_space< IntType > * evts
                            , device_allele_space< RealType, IntType, unordered_tag > * out_space ) {

}

#endif  // DEVICE_ALLELE_SPACE_UNIT_ORDERED_KERNELS_HPP_
