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
#ifndef SIMPLE_MUTATION_GENERATOR_HPP_
#define SIMPLE_MUTATION_GENERATOR_HPP_

#include <cuda.h>

//#include "clotho/cuda/data_spaces/allele_space/device_allele_space.hpp"
#include "clotho/cuda/data_spaces/event_space/device_event_space.hpp"
#include "clotho/cuda/distributions/poisson_distribution.hpp"

/*
template < class StateType, class RealType, class IntType >
__global__ void _simple_mutation_generator( StateType * states
                                            , device_event_space< IntType, unordered_tag > * mut_events
                                            , poisson_cdf< RealType, 32 > * dist
                                            , unsigned int sequence_size ) {
    typedef device_event_space< IntType, unordered_tag > space_type;
    typedef typename space_type::pointer    event_ptr;

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    __shared__ RealType cdf[ 32 ];

    unsigned int max_k = dist->max_k;
    cdf[ tid ] = dist->_cdf[ tid ];
    __syncthreads();

    StateType local_state = states[ tid ];

    event_ptr evt = mut_events->event_count;
    __syncthreads();

    unsigned int i = tid;
    unsigned int count = 0;
    while ( i < sequence_size ) {
        RealType x = curand_uniform( &local_state );
        unsigned int _c = _find_poisson_maxk32( cdf, x, max_k );

        evt[ i ] = _c;

        count += _c;
        i += (blockDim.x * blockDim.y);
    }
    __syncthreads();

    for( i = 1; i < 32; i <<= 1 ) {
        unsigned int _c = __shfl_up( count, i );
        count += ( tid >= i ) * _c;
    }

    if( tid == 31 ) {
        // has total number of mutations generated
        mut_events->total = count;    
    }

    states[ tid ] = local_state;
}

template < class StateType, class RealType, class IntType >
__global__ void _simple_mutation_generator( StateType * states
                                            , device_event_space< IntType, unit_ordered_tag< IntType > > * mut_events
                                            , poisson_cdf< RealType, 32 > * dist
                                            , unsigned int sequence_size ) {

    typedef device_event_space< IntType, unit_ordered_tag< IntType > > space_type;
    typedef typename space_type::pointer    event_ptr;

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    __shared__ RealType cdf[ 32 ];

    unsigned int max_k = dist->max_k;
    cdf[ tid ] = dist->_cdf[ tid ];
    __syncthreads();

    StateType local_state = states[ tid ];

    event_ptr evt = mut_events->event_count;
    __syncthreads();

    unsigned int i = tid, count = 0;
    while ( sequence_size-- ) {
        RealType x = curand_uniform( &local_state );
        unsigned int _c = _find_poisson_maxk32( cdf, x, max_k );

        evt[ i ] = _c;

        count += _c;
        i += (blockDim.x * blockDim.y);
    }
    __syncthreads();

    mut_events->bin_summary[tid] = count;

    for( i = 1; i < 32; i <<= 1 ) {
        unsigned int _c = __shfl_up( count, i );
        count += ( tid >= i ) * _c;
    }

    if( tid == 31 ) {
        // has total number of mutations generated
        mut_events->total = count;    
    }

    states[ tid ] = local_state;
}

template < class StateType, class RealType, class IntType >
__global__ void _simple_mutation_generator2( StateType * states
                                            , device_event_space< IntType, unordered_tag > * mut_events
                                            , poisson_cdf< RealType, 32 > * dist
                                            , unsigned int sequence_size ) {
    typedef device_event_space< IntType, unordered_tag > space_type;
    typedef typename space_type::pointer    event_ptr;

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    __shared__ RealType cdf[ 32 ];

    unsigned int max_k = dist->max_k;
    cdf[ tid ] = dist->_cdf[ tid ];
    __syncthreads();

    StateType local_state = states[ tid ];

    event_ptr evt = mut_events->event_count;
    __syncthreads();

    unsigned int i = tid;
    unsigned int count = 0;

    unsigned int padded_size = ((sequence_size & 31) ? ((sequence_size & (~31)) + 32) : sequence_size);  // ==> % 32

    while ( i < padded_size ) {
        RealType x = curand_uniform( &local_state );
        unsigned int _c = _find_poisson_maxk32( cdf, x, max_k );

        _c *= (i < sequence_size);

        for( unsigned int j = 1; j < 32; j <<= 1 ) {
            unsigned int t = __shfl_up( _c, j );
            _c += (tid >= j) * t;
        }

        count += _c;

        if( i < sequence_size ) {
            evt[ i ] = count;
        }
        __syncthreads();
        i += (blockDim.x * blockDim.y);

        unsigned int t = __shfl(count, 31);
        count = t;
    }
    __syncthreads();

    if( tid == 31 ) {
        // has total number of mutations generated
        mut_events->total = count;
    }

    states[ tid ] = local_state;
}*/

template < class StateType, class RealType, class IntType >
__global__ void make_event_distribution_kernel( StateType * states
                                            , device_event_space< IntType, unordered_tag > * _events
                                            , RealType lambda
                                            , unsigned int N ) {
    typedef StateType state_type;
    typedef RealType real_type;
 
    if( (threadIdx.y * blockDim.x + threadIdx.x) == 0 ) {
        lambda *= (real_type) N;
        unsigned int n = curand_poisson( states, lambda );

        _events->total = n;
    }
}

template < class StateType, class RealType, class IntType >
__global__ void make_event_distribution_kernel( StateType * states
                                                , device_event_space< IntType, unit_ordered_tag< IntType > > * _events
                                                , RealType lambda
                                                , unsigned int N ) {
    typedef StateType   state_type;
    typedef RealType    real_type;
    typedef unit_ordered_tag< IntType > ordered_tag_type;

    if( blockIdx.y * gridDim.x + blockIdx.x != 0 ) return;

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int tpb = blockDim.x * blockDim.y;

    assert( tpb == 32 );
    assert( ordered_tag_type::OBJECTS_PER_UNIT == 32 );

    lambda /= ((real_type) ordered_tag_type::OBJECTS_PER_UNIT );
    lambda *= ((real_type) N );

    state_type local_state = states[ tid ];

    unsigned int count = curand_poisson( &local_state, lambda );
    __syncthreads();

    _events->bin_summary[ tid ] = count;
    __syncthreads();

    states[ tid ] = local_state;

    for( unsigned int i = 1; i < 32; i <<= 1 ) {
        unsigned int c = __shfl_up( count, i );
        count += ((tid >= i ) * c );
    }

    if( tid == 31 ) {
        _events->total = count;
    }
}


#endif  // SIMPLE_MUTATION_GENERATOR_HPP_
