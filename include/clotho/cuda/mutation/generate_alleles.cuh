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
#ifndef GENERATE_ALLELE_CUH_
#define GENERATE_ALLELE_CUH_

#include "clotho/cuda/allele_space/device_allele_space.hpp"
#include "clotho/cuda/mutation/mutation_space.hpp"

template < class RealType >
__device__ unsigned int get_bin_index( RealType x, unsigned int bin_count = 32 ) {
    unsigned int idx = (unsigned int) floor( x * (RealType) bin_count );
    return idx;
}

template < class StateType, class RealType, class IntType >
__global__ void generate_alleles( StateType * states
                                , device_allele_space< RealType, IntType, partially_ordered_tag > * aspace
                                , mutation_space< RealType, IntType > * mspace ) {
    IntType tid = threadIdx.y * blockDim.x + threadIdx.x;

    StateType local_state = states[ blockIdx.x * blockDim.x * blockDim.y + tid ];

    __shared__ IntType free_bits[ sizeof(IntType) * 8 ];

    // generate a random location from a uniform distribution
    while( N ) {
        RealType loc = curand_uniform( &local_state );
        unsigned int bin_idx = get_bin_index( loc );

        

        N -= (( N > blockDim.x*blockDim.y ) ? (blockDim.x * blockDim.y) : N);
    }


    states[ blockIdx.x * blockDim.x * blockDim.y + tid ] = local_state;
}

#endif  // GENERATE_ALLELE_CUH_
