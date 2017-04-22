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
    memcpy( out_space->locations, in_space->locations, in_space->capacity * sizeof( typename device_allele_space< RealType >::real_type ) );
    memcpy( out_space->neutral, in_space->neutral, in_space->capacity * sizeof( typename device_allele_space< RealType >::real_type ) );
}

template < class RealType >
__device__ void _update_space( device_weighted_allele_space< RealType > * in_space
                            , device_weighted_allele_space< RealType > * out_space ) {
    _update_space( (device_allele_space< RealType > * ) in_space, (device_allele_space< RealType > * ) out_space );

    memcpy( out_space->weights, in_space->weights, in_space->capacity * sizeof( typename device_allele_space< RealType >::real_type ) );
}

template < class AlleleSpaceType, class IntType >
__global__ void _merge_space( AlleleSpaceType * in_alleles
                            , device_free_space< IntType, unordered_tag > * in_free
                            , device_event_space< IntType, unordered_tag > * new_muts
                            , device_free_space< IntType, unordered_tag > * out_free
                            , AlleleSpaceType * out_alleles ) {

    unsigned int N = in_free->total;    // total number of free elements from parent allele space
    unsigned int M = new_muts->total;   // total number of new mutations to be added to offspring population

    unsigned int P = in_alleles->capacity; // current allocated size of parent population

    unsigned int C = P + (( N < M ) ? (M - N) : 0);

//    printf( "Resize Population: %d = %d + |min(0, (%d - %d))|\n", C, P, N, M );

    _resize_space_impl(out_alleles, C );  // resize offspring allele space relative to parent space

    C = out_alleles->capacity;  // the size of the offspring allele space
//    printf( "Padded Population: %d\n", C );

    _resize_space_impl( out_free, C );     // resize offspring free space relative to offsrping allele space
 
    if( P > 0 ) {
        //printf( "Updating spaces \n" );
        _update_space( in_alleles, out_alleles );
        _update_space( in_free, out_free );
    }
}

template < class AlleleSpaceType, class IntType >
__global__ void expand_spaces_kernel( AlleleSpaceType * in_alleles
                            , device_free_space< IntType, unordered_tag > * in_free
                            , device_event_space< IntType, unordered_tag > * new_muts
                            , device_free_space< IntType, unordered_tag > * out_free
                            , AlleleSpaceType * out_alleles ) {

    unsigned int N = in_free->total;    // total number of free elements from parent allele space
    unsigned int M = new_muts->total;   // total number of new mutations to be added to offspring population

    unsigned int P = in_alleles->capacity; // current allocated size of parent population

    unsigned int C = P + (( N < M ) ? (M - N) : 0);

    _resize_space_impl(out_alleles, C );  // resize offspring allele space relative to parent space

    C = out_alleles->capacity;  // the size of the offspring allele space

    _resize_space_impl( out_free, C );     // resize offspring free space relative to offsrping allele space
}

template < class RealType >
__global__ void update_space_kernel( device_allele_space< RealType > * in_alleles
                                     , device_allele_space< RealType > * out_alleles ) {

    if( blockIdx.y * gridDim.x + blockIdx.x > 0 ) return;

    unsigned int N = in_alleles->capacity;
    unsigned int M = out_alleles->capacity;

    assert( N <= M );

    unsigned int tpb = blockDim.x * blockDim.y;
    unsigned int idx = threadIdx.y * blockDim.x + threadIdx.x;

    RealType * in_data = in_alleles->location;
    RealType * out_data = out_alleles->location;

    while( idx < N ) {
        out_data[idx] = in_data[idx];
        idx += tpb;
    }
    __syncthreads();

    while( idx < M ) {
        out_data[idx] = (RealType)0;
        idx += tpb;
    }
}

template < class RealType >
__global__ void update_space_kernel( device_weighted_allele_space< RealType > * in_alleles
                                     , device_weighted_allele_space< RealType > * out_alleles ) {

    
    unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;
    if( bid > 2 ) return;

    unsigned int N = in_alleles->capacity;
    unsigned int M = out_alleles->capacity;

    assert( N <= M );

    unsigned int tpb = blockDim.x * blockDim.y;
    unsigned int idx = threadIdx.y * blockDim.x + threadIdx.x;

    RealType * in_data, * out_data;

    if( bid == 0 ) {
        in_data = in_alleles->locations;
        out_data = out_alleles->locations;
    } else if( bid == 1 ) {
        in_data = in_alleles->neutral;
        out_data = out_alleles->neutral;
    } else {
        in_data = in_alleles->weights;
        out_data = out_alleles->weights;
    }

    // Copy in space to out space
    while( idx < N ) {
        out_data[idx] = in_data[idx];
        idx += tpb;
    }
    __syncthreads();

    // clear out space that is unused
    while( idx < M ) {
        out_data[idx] = (RealType)0;
        idx += tpb;
    }
}

#endif  // DEVICE_ALLELE_SPACE_UNORDERED_KERNELS_HPP_
