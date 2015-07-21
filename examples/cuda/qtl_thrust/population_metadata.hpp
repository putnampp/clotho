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
#ifndef POPULATION_METADATA_HPP_
#define POPULATION_METADATA_HPP_

#include <cuda.h>

/*
struct fixed_op {
    inline __device__ unsigned int operator()( unsigned int a, unsigned int b ) {
        return (a & b);
    }
};

struct variable_op {
    inline __device__ unsigned int operator()( unsigned int a, unsigned int b ) {
        return (a | b);
    }
};

template < class OP >
struct population_metadata {

    __global__ void eval( unsigned int * pop
                        , unsigned int rows
                        , unsigned int cols
                        , unsigned int * metadata ) {

        unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

        __shared__ unsigned int sPop[ MAX_THREADS ];
        __shared__ unsigned int sMeta[ MAX_THREADS ];

        unsigned int s_offset = ((threadIdx.y + 1) * blockDim.x + threadIdx.x)

        OP oper;

        unsigned int i = 0;
        while( i < rows ) {
            sPop[ tid ] = pop[ offset ];
            __syncthreads();

            unsigned int d = sPop[tid], e = sPop[s_offset];
            unsigned int res = sMeta[tid];
            __syncthreads();

            d = oper.eval( d, e);
            res = oper.eval( res, d);
            __syncthreads();

            sMeta[tid] = res;
            __syncthreads();

            i += blockDim.y;
            offset += blockDim.y * cols;
        }
    }
};*/

__global__ void update_population_metadata( unsigned int * pop
                                    , unsigned int rows
                                    , unsigned int cols
                                    , unsigned int * free
                                    , unsigned int * lost
                                    , unsigned int * fixed );

__global__ void update_population_fixed( unsigned int * pop
                                    , unsigned int rows
                                    , unsigned int cols
                                    , unsigned int * fixed );

__global__ void update_population_lost( unsigned int * pop
                                    , unsigned int rows
                                    , unsigned int cols
                                    , unsigned int * fixed );

__global__ void update_population_free( unsigned int * fixed
                                    , unsigned int * lost
                                    , unsigned int * free
                                    , unsigned int cols
                                    );
#endif  // POPULATION_METADATA_HPP_
