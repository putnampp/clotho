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
#ifndef CUDA_BERNOULLI_DISTRIBUTION_HPP_
#define CUDA_BERNOULLI_DISTRIBUTION_HPP_

#include <cuda.h>

template < class StateType, class SpaceType, class RealType >
__global__ void bernoulli_kernel( StateType * states, SpaceType * space, RealType bias_rate ) {
    assert(0 <= bias_rate && bias_rate <= 1.0 );
    if( bias_rate == 0.) return; // do not change

    typedef typename SpaceType::value_type  int_type;
    unsigned int N = space->size;
    if( N == 0 ) return;

    int_type * d = space->data;

    unsigned int idx = threadIdx.y * blockDim.x + threadIdx.x;
    StateType local_state = states[ idx ];

    unsigned int tpb = (blockDim.x * blockDim.y);
    unsigned int pad_N = N / tpb;
    pad_N += ((N % tpb) ? 1 : 0);
    pad_N *= tpb;

    while( idx < pad_N ) {
        RealType x = curand_uniform( &local_state );

        if( idx < N && (x >= bias_rate) ) {
            d[ idx ] = 1;
        }
        __syncthreads();

        idx += tpb;
    }

    idx = threadIdx.y * blockDim.x + threadIdx.x;
    states[ idx ] = local_state;
}

template < class StateType, class SpaceType, class RealType, class ShiftType >
__global__ void inline_bernoulli_linear_shift_kernel( StateType * states, SpaceType * space, RealType bias_rate, ShiftType shift ) {
    assert(0 <= bias_rate && bias_rate <= 1.0 );
    if( bias_rate == 0.) return; // do not change

    typedef typename SpaceType::value_type  int_type;
    unsigned int N = space->size;
    if( N == 0 ) return;

    int_type * d = space->data;

    unsigned int idx = threadIdx.y * blockDim.x + threadIdx.x;
    StateType local_state = states[ idx ];

    unsigned int tpb = (blockDim.x * blockDim.y);
    unsigned int pad_N = N / tpb;
    pad_N += ((N % tpb) ? 1 : 0);
    pad_N *= tpb;

    while( idx < pad_N ) {
        RealType x = curand_uniform( &local_state );

        if( idx < N && (x >= bias_rate) ) {
            d[ idx ] += shift;
        }
        __syncthreads();

        idx += tpb;
    }

    idx = threadIdx.y * blockDim.x + threadIdx.x;
    states[ idx ] = local_state;
}
#endif  // CUDA_BERNOULLI_DISTRIBUTION_HPP_
