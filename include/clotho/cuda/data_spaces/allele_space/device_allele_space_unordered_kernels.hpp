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

/*
template < class RealType, class IntType >
__global__ void _update_free_count( device_allele_space< RealType, IntType, unordered_tag > * aspace ) {
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    IntType * flist = aspace->free_list;
    IntType N = aspace->free_list_size();
    popcountGPU< IntType > pc;

    unsigned int _count = 0;
    while( tid < N ) {
        IntType c = flist[tid];

        _count += pc.evalGPU( c );
        tid += 32;
    }
    __syncthreads();

    for( unsigned int i = 1; i < 32; i <<= 1 ) {
        IntType t = __shfl_up( _count, i );
        _count += (( (tid & 31) >= i ) * t);
    }

    if( (tid & 31) == 31) {
        aspace->free_count = _count;
    }
}*/
/*
template < class RealType, class IntType >
__global__ void _merge_allele_space( device_allele_space< RealType, IntType, unordered_tag > * aspace
                                    , device_allele_space< RealType, IntType, unordered_tag > * bspace
                                    , device_allele_space< RealType, IntType, unordered_tag > * output ) {
    typedef device_allele_space< RealType, IntType, unordered_tag > space_type;

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    space_type A = *aspace;
    space_type B = *bspace;
    space_type O = *output;

    unsigned int min_out_size = (A.size - A.free_count);
    min_out_size += (B.size - B.free_count);

    if( (tid == 0) && (O.size < min_out_size) ) {
        // resize output
        _resize_allele_space( &O, min_out_size );
    }
    __syncthreads();

    // persist output allele space
    *output = O;
}*/

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
//        memcpy( out_space->free_list, local.free_list, local.free_list_size() * sizeof( typename device_allele_space< RealType, IntType, unordered_tag >::int_type ));
    }
}
#endif  // DEVICE_ALLELE_SPACE_UNORDERED_KERNELS_HPP_
