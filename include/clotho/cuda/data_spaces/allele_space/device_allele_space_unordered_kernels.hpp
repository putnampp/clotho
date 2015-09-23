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

template < class RealType >
__device__ void _update_space( device_allele_space< RealType > * in_space
                            , device_allele_space< RealType > * out_space ) {
    memcpy( out_space->locations, in_space->locations, in_space->size * sizeof( typename device_allele_space< RealType >::real_type ) );
}

template < class RealType >
__device__ void _update_space( device_weighted_allele_space< RealType > * in_space
                            , device_weighted_allele_space< RealType > * out_space ) {
    _update_space( (device_allele_space< RealType > * ) in_space, (device_allele_space< RealType > * ) out_space );

    memcpy( out_space->weights, in_space->weights, in_space->size * sizeof( typename device_allele_space< RealType >::real_type ) );
}

/*
template < class RealType, class IntType >
__global__ void _merge_space( device_allele_space< RealType > * in_space
                            , device_free_space< IntType, unordered_tag > * fspace
                            , device_event_space< IntType, unordered_tag > * evts
                            , device_allele_space< RealType > * out_space ) {

    typedef device_allele_space< RealType > space_type;
    space_type local_in = *in_space;

    unsigned int N = local_in.size - fspace->total;
    N += evts->total;

    if( local_in.size > N ) { N = local_in.size; }

    _resize_space_impl( out_space, N ); // will pad space based upon ALIGNMENT_SIZE

    if( local_in.locations ) {
        memcpy( out_space->locations, local_in.locations, local_in.size * sizeof( typename space_type::real_type ) );
    }
}
*/

/*
template < class RealType, class IntType >
__global__ void _merge_space( device_weighted_allele_space< RealType > * in_space
                            , device_free_space< IntType, unordered_tag > * fspace
                            , device_event_space< IntType, unordered_tag > * evts
                            , device_free_space< IntType, unordered_tag > * ofspace
                            , device_allele_space< RealType > * out_space ) {

    typedef device_allele_space< RealType > space_type;
    space_type local_in = *in_space;

    unsigned int N = local_in.size - fspace->total;
    N += evts->total;

    if( local_in.size > N ) {
        N = local_in.size; 
    } else {
        _resize_space_impl( out_space, N ); // will pad space based upon ALIGNMENT_SIZE
        _resize_space_impl( out_space, ofspace );
    }

//    memcpy( out_space->locations, local_in.locations, local_in.size * sizeof( typename space_type::real_type ) );
//    memcpy( out_space->weights, local_in.weights, local_in.size * sizeof( typename space_type::real_type ) );

    _update_space( in_space, out_space );
    _update_space( fspace, ofspace );    
}*/

template < class RealType, class IntType >
__global__ void _merge_space( device_weighted_allele_space< RealType > * in_alleles
                            , device_free_space< IntType, unordered_tag > * in_free
                            , device_event_space< IntType, unordered_tag > * new_muts
                            , device_free_space< IntType, unordered_tag > * out_free
                            , device_weighted_allele_space< RealType > * out_alleles ) {

    unsigned int N = in_free->total;    // total number of free elements from parent allele space
    unsigned int M = new_muts->total;   // total number of new mutations to be added to offspring population

    unsigned int P = in_alleles->size; // current allocated size of parent population

    unsigned int C = P + (( N < M ) ? (M - N) : 0);

//    printf( "Resize Population: %d = %d + |min(0, (%d - %d))|\n", C, P, N, M );

    _resize_space_impl(out_alleles, C );  // resize offspring allele space relative to parent space

    C = out_alleles->size;  // the size of the offspring allele space
//    printf( "Padded Population: %d\n", C );

    _resize_space_impl( out_free, C );     // resize offspring free space relative to offsrping allele space
    
    _update_space( in_alleles, out_alleles );
    _update_space( in_free, out_free );
}

#endif  // DEVICE_ALLELE_SPACE_UNORDERED_KERNELS_HPP_
