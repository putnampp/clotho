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
#ifndef DEVICE_ALLELE_SPACE_UNORDERED_KERNELS_HPP_
#define DEVICE_ALLELE_SPACE_UNORDERED_KERNELS_HPP_

#include "clotho/cuda/data_spaces/allele_space/device_allele_space_def.hpp"
#include "clotho/cuda/data_spaces/allele_space/device_allele_space_kernel_api.hpp"

#include "clotho/cuda/data_spaces/tags/unordered_tag.hpp"
#include "clotho/cuda/popcount_kernel.h"

#include <stdlib.h>

template < class RealType, class IntType >
__global__ void _merge_space( device_allele_space< RealType/*, IntType, unordered_tag*/ > * in_space
                            , device_free_space< IntType, unordered_tag > * fspace
                            , device_event_space< IntType, unordered_tag > * evts
                            , device_allele_space< RealType/*, IntType, unordered_tag*/ > * out_space ) {

    typedef device_allele_space< RealType/*, IntType, unordered_tag*/ > space_type;
    space_type local = *in_space;

    unsigned int N = local.size - fspace->total;
    N += evts->total;

    if( local.size > N ) { N = local.size; }

    _resize_space_impl( out_space, N ); // will pad space based upon ALIGNMENT_SIZE

    if( local.locations ) {
        memcpy( out_space->locations, local.locations, local.size * sizeof( typename space_type::real_type ) );
    }
}
#endif  // DEVICE_ALLELE_SPACE_UNORDERED_KERNELS_HPP_
