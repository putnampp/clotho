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
#ifndef GENERATE_MUTATION_KERNEL_IMPL_HPP_
#define GENERATE_MUTATION_KERNEL_IMPL_HPP_

#include "clotho/cuda/data_spaces/allele_space/device_allele_space.hpp"
#include "clotho/cuda/data_spaces/event_space/device_event_space.hpp"
#include "clotho/cuda/data_spaces/free_space/device_free_space.hpp"

template < class StateType, class AlleleSpaceType, class IntType, class OrderTag >
__global__ void _generate_mutation_kernel( StateType * states
                                        , device_free_space< IntType, OrderTag > * fspace
                                        , device_event_space< IntType, OrderTag > * events
                                        , AlleleSpaceType * alleles ) {
    typedef StateType                               state_type;

    typedef device_free_space< IntType, OrderTag >  free_space_type;
    typedef device_event_space< IntType, OrderTag > event_space_type;

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    unsigned int N = events->total;
    unsigned int i = tid;
    unsigned int * fmap = fspace->free_map;

    state_type local_state = states[tid];

    OrderTag * otag = NULL;

    while( i < N ) {
        unsigned int idx = fmap[ i ];

        _generate_random_allele( &local_state, alleles, idx, otag );
        
        i += (blockDim.x * blockDim.y);
    }

    states[tid] = local_state;
}

/*
template < class StateType, class RealType, class IntType, class OrderTag >
__global__ void _generate_mutation_kernel( StateType * states
                                        , device_free_space< IntType, OrderTag > * fspace
                                        , device_event_space< IntType, OrderTag > * events
                                        , device_weighted_allele_space< RealType > * alleles ) {
    typedef StateType                               state_type;

    typedef device_free_space< IntType, OrderTag >  free_space_type;

    typedef device_event_space< IntType, OrderTag > event_space_type;

    typedef device_weighted_allele_space< RealType >  allele_space_type;
    typedef typename allele_space_type::real_type     real_type;

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    unsigned int N = events->total;

    unsigned int i = tid;

//    allele_space_type local_space = * alleles;

    real_type * locs = alleles->locations;
    real_type * wghts = alleles->weights;
    
    unsigned int * fmap = fspace->free_map;

    state_type local_state = states[tid];

    while( i < N ) {
        unsigned int idx = fmap[ i ];

        real_type x = curand_uniform( &local_state );
        real_type y = curand_normal( &local_state ); 

        locs[idx] = x;
        wghts[idx] = y;
        
        i += (blockDim.x * blockDim.y);
    }

    states[tid] = local_state;
}*/

#endif  // GENERATE_MUTATION_KERNEL_IMPL_HPP_

