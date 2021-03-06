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
#ifndef CURAND_HELPER_HPP_
#define CURAND_HELPER_HPP_

#include <string>
#include <curand.h>
#include <curand_kernel.h>
#include <curand_mtgp32.h>

namespace clotho {
namespace cuda {

template < class StateType >
__global__ void setup_state_kernel( StateType * states, unsigned long long seed ) {
    int id = blockIdx.x * blockDim.x * blockDim.y + threadIdx.y * blockDim.x + threadIdx.x;

    curand_init( seed, id, 0, &states[id] );
}

template < class State >
struct curand_helper {
    typedef State               state_type;
    typedef unsigned long long  seed_type;

    static const std::string StateName;

    static void make_states( state_type *& states, seed_type seed, unsigned int blocks, unsigned int threads );

    static void cleanup_states( state_type * states );
};

template <>
const std::string curand_helper< curandStateMtgp32_t >::StateName = "MTGP32";

template <>
void curand_helper< curandStateMtgp32_t >::make_states( state_type *& states, seed_type seed, unsigned int blocks, unsigned int threads ) {
    assert(false);
}

template <>
void curand_helper< curandStateMtgp32_t >::cleanup_states( state_type * states ) {
    cudaFree( states );
}

template <>
const std::string curand_helper< curandStateXORWOW >::StateName = "XORWOW";

template <>
void curand_helper< curandStateXORWOW >::cleanup_states( state_type * states ) {
    if( states != NULL )
        cudaFree( states );
}

template <>
void curand_helper< curandStateXORWOW >::make_states( state_type *& states, seed_type seed, unsigned int blocks, unsigned int threads ) {
//    if( states != NULL ) {
//        cleanup_states( states );
//    }

    assert( cudaMalloc( (void **) &states, blocks * threads * sizeof(state_type) ) == cudaSuccess );

    setup_state_kernel<<< blocks, threads >>>( states, seed );
}

}   // namespace cuda
}   // namespace clotho


#endif  // CURAND_HELPER_HPP_
