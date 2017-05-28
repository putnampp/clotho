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
#ifndef SEQUENCE_WEIGHT_KERNEL_HPP_
#define SEQUENCE_WEIGHT_KERNEL_HPP_

template < class IntType >
__global__ void evaluate_sequence_weights( IntType * pop, unsigned int * weight, unsigned int width ) {
    assert( blockDim.x == 32 );

    unsigned int b_idx = threadIdx.y;
    
    IntType mask = (1 << threadIdx.x);

    __shared__ unsigned int sBuffer[ 32 ];

    if( threadIdx.y == 0 ) {
        sBuffer[ threadIdx.x ] = 0;
    }
    __syncthreads();

    unsigned int N = 0;
    while( b_idx < width ) {
        IntType b = pop[ blockIdx.x * width + b_idx ];
        if( b & mask ) {
            N += 1;
        }
        b_idx += blockDim.y;
    }
    __syncthreads();

    for( unsigned int i = 1; i < 32; i <<= 1) {
        unsigned int _N = __shfl_up( N, i );
        N += (( threadIdx.x > i) ? _N : 0);
    }

    if( threadIdx.x == 31 ) {
        sBuffer[ threadIdx.y ] = N;
    }
    __syncthreads();

    N = sBuffer[ threadIdx.x ];
    for( unsigned int i = 1; i < 32; i <<= 1 ) {
        unsigned int _N = __shfl_up( N, i );
        N += (( threadIdx.x > i) ? _N : 0 );
    }

    if( threadIdx.y == 0 && threadIdx.x == 31 ) {
        weight[ blockIdx.x ] = N;
    }
}

#endif  // SEQUENCE_WEIGHT_KERNEL_HPP_
