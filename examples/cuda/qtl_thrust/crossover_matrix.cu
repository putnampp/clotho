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
#include "crossover_matrix.hpp"

const unsigned int MAX_THREADS = 1024;

__global__ void generate_crossover_matrix( double * rand_pool, unsigned int * event_list, double * alleles, unsigned int * sequences, dim3 max_dims ) {
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
//    unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;
//    unsigned int tcount = blockDim.x * blockDim.y * blockDim.z;

    // all threads within a block will operate upon a contiguous set of events and alleles
    __shared__ double lEvents[ MAX_THREADS ];       //  8K
    __shared__ double lAlleles[ MAX_THREADS ];      //  8K
    __shared__ unsigned int lMask[ MAX_THREADS ];   //  4K
    __shared__ unsigned int lCounts[ MAX_THREADS ];  //  4K 
                                                    // 24K shared memory

    // copy the alleles that this block is going to be inspecting
    unsigned int aidx = blockIdx.x * MAX_THREADS + threadIdx.x;
    lAlleles[ tid ] = ((aidx < max_dims.x ) ? alleles[aidx] : -1.0);

    __syncthreads();

    unsigned int seq_per_block = max_dims.y / gridDim.y;

    unsigned int seq_start = gridDim.y * seq_per_block, seq_end = seq_start + seq_per_block;

    if( seq_end > max_dims.y ) seq_end = max_dims.y;

    unsigned int c = 0;
    unsigned int mask = (1 << threadIdx.x);
    while( seq_start < seq_end ) {
        if( c == 0 ) {
            // refresh event count list
            c = seq_end - seq_start;
            if( c > MAX_THREADS ) { c = MAX_THREADS; }

            lCounts[ tid ] = (( tid < c ) ? event_list[ seq_start + tid ] : 0);
            --c;
        }
        __syncthreads();

        lMask[tid] = 0; // reset masks for sequence
        __syncthreads();

        // number of events for sequence
        unsigned int eEnd = lCounts[ c ], eStart = lCounts[ --c ];
        unsigned int N = 0;
        while( eStart < eEnd ) {
            if( N == 0 ) {
                // refresh event list

                // current method performs a 'lookahead copy'.
                // general question: is it better to perform a block copy from global memory to local memory,
                // or to simply access global memory as necessary?

                // Basic concern is that when N is small there will be a lot of sleeping threads
                N = eEnd - eStart;
                if( N > MAX_THREADS ) { N = MAX_THREADS; }
                lEvents[ tid ] = ((tid < N ) ? rand_pool[ eStart + tid ] : 0);
            }
            __syncthreads();

            eStart += N;

            while( N-- ) {
                lMask[tid] ^= (( lAlleles[tid] >= 0 && lAlleles[tid] < lEvents[N]) ? mask : 0 );
            }
            __syncthreads();

            N = 0;
        }

        // collapse masks for sequence
        if( threadIdx.x & 1 == 0 ) {
            // half warp
            // merge odds and evens: 0 < 1; 2 < 3; 4 < 5 ...
            lMask[tid] |= lMask[tid + 1];
        }

        __syncthreads();

        if( threadIdx.x & 3 == 0 ) {
            // quarter warp
            // merge 0 < 2;4 < 6;8 < 10...
            lMask[tid] |= lMask[tid + 2];
        }

        __syncthreads();

        if( threadIdx.x & 7 == 0 ) {
            // eighth warp
            // merge 0 < 4; 8 < 12; 16 < 20; 24 < 28
            lMask[tid] |= lMask[tid + 4];
        }

        __syncthreads();

        if( threadIdx.x & 15 == 0 ) {
            // sixth warp
            // merge 0 < 8; 16 < 24
            lMask[tid] |= lMask[tid + 8];
        }

        __syncthreads();

        // copy sequence to output
        if( threadIdx.x & 31 == 0 ) {
            // thirty-second warp (1-thread)
            // merge 0 < 16
            lMask[ tid ] |= lMask[ tid + 16 ];
        }

        __syncthreads();

        if( tid < 32 ) {
            // use a single warp to copy between shared and global memory
            sequences[ seq_start * max_dims.x / 32 + tid ] = lMask[ tid * 32 ];
        }
        
        __syncthreads();
        ++seq_start;
    }
}

__global__ void generate_crossover_matrix2( double * rand_pool, unsigned int * event_list, double * alleles, unsigned int * sequences, dim3 max_dims ) {
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
//    unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;
//
//    unsigned int tcount = blockDim.x * blockDim.y * blockDim.z;

    // all threads within a block will operate upon a contiguous set of events and alleles
    __shared__ double lEvent;
    __shared__ double lAlleles[ MAX_THREADS ];      //  8K
    __shared__ unsigned int lMask[ MAX_THREADS ];   //  4K
    __shared__ unsigned int lCounts[ MAX_THREADS ];  //  4K 
                                                    // 24K shared memory

    // copy the alleles that this block is going to be inspecting
    unsigned int aidx = blockIdx.x * MAX_THREADS + threadIdx.x;
    lAlleles[ tid ] = ((aidx < max_dims.x ) ? alleles[aidx] : -1.0);

    __syncthreads();

    unsigned int seq_per_block = max_dims.y / gridDim.y;

    unsigned int seq_start = gridDim.y * seq_per_block, seq_end = seq_start + seq_per_block;

    if( seq_end > max_dims.y ) seq_end = max_dims.y;

    unsigned int c = 0;
    unsigned int mask = (1 << threadIdx.x);
    while( seq_start < seq_end ) {
        if( c == 0 ) {
            // refresh event count list
            c = seq_end - seq_start;
            if( c > MAX_THREADS ) { c = MAX_THREADS; }

            lCounts[ tid ] = (( tid < c ) ? event_list[ seq_start + tid ] : 0);
            --c;
        }
        __syncthreads();

        lMask[tid] = 0; // reset masks for sequence
        __syncthreads();

        // number of events for sequence
        unsigned int eEnd = lCounts[ c ], eStart = lCounts[ --c ];
        while( eStart < eEnd ) {
            if( tid == 0 ) {
                // use one thread to
                // get the next event from GLOBAL memory
                lEvent = rand_pool[ eStart ];
            }
            __syncthreads();

            lMask[tid] ^= (( lAlleles[tid] >= 0 && lAlleles[tid] < lEvent) ? mask : 0 );
            __syncthreads();
            ++eStart;
        }

        // collapse masks for sequence
        if( threadIdx.x & 1 == 0 ) {
            // half warp
            // merge odds and evens: 0 < 1; 2 < 3; 4 < 5 ...
            lMask[tid] |= lMask[tid + 1];
        }

        __syncthreads();

        if( threadIdx.x & 3 == 0 ) {
            // quarter warp
            // merge 0 < 2;4 < 6;8 < 10...
            lMask[tid] |= lMask[tid + 2];
        }

        __syncthreads();

        if( threadIdx.x & 7 == 0 ) {
            // eighth warp
            // merge 0 < 4; 8 < 12; 16 < 20; 24 < 28
            lMask[tid] |= lMask[tid + 4];
        }

        __syncthreads();

        if( threadIdx.x & 15 == 0 ) {
            // sixth warp
            // merge 0 < 8; 16 < 24
            lMask[tid] |= lMask[tid + 8];
        }

        __syncthreads();

        // copy sequence to output
        if( threadIdx.x & 31 == 0 ) {
            // thirty-second warp (1-thread)
            // merge 0 < 16
            lMask[ tid ] |= lMask[ tid + 16 ];
        }

        __syncthreads();

        if( tid < 32 ) {
            // use a single warp to copy between shared and global memory
            sequences[ seq_start * max_dims.x + tid ] = lMask[ tid * 32 ];
        }
        
        __syncthreads();
        ++seq_start;
    }
}
