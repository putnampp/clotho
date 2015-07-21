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
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either expcount_ptrs or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License.
#include "scatter_mutations.hpp"

#include "popcount_kernel.h"

const unsigned int BLOCK_PER_STRIDE = 32;
const unsigned int MAX_THREADS = 1024;

inline __device__ void scan_up( unsigned int val, unsigned int * d ) {
#pragma unroll
    for(unsigned int i = 1; i < BLOCK_PER_STRIDE; i <<= 1 ) {
        unsigned int n = __shfl_up(val, i, BLOCK_PER_STRIDE );
        if( threadIdx.x >= i ) val += n;
    }

    if( threadIdx.y == 0 ) {
        d[ threadIdx.x + 1 ] = val;
    }
    __syncthreads();
}

__global__ void scatter_mutations( double * rands
                                 , double * alleles
                                 , unsigned int * free_list
                                 , unsigned int * sequences
                                 , unsigned int * mut_events
                                 , unsigned int seq_count
                                 , unsigned int allele_count
                                 , unsigned int * dbg
                                ) {

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

//    __shared__ double lMutations[ MAX_THREADS ];                // 8 * 1024
//    __shared__ double lAllele[ MAX_THREADS ];                   // 8 * 1024
    __shared__ unsigned int lFree[ BLOCK_PER_STRIDE ];         // 4 * 32
    __shared__ unsigned int lCounts[ BLOCK_PER_STRIDE + 1 ];   // 4 * (32 + 1)
//    __shared__ unsigned int lEvents[ MAX_THREADS ];             // 4 * 1024
//    __shared__ unsigned int lSeq[ BLOCK_PER_STRIDE ];           // 4 * 32
//    __shared__ unsigned int exp_free[ MAX_THREADS ];
                                                                // < 24 K
    popcountGPU< unsigned int > pc;

    // make sure lCounts is cleared
    if( tid <= BLOCK_PER_STRIDE ) {
        lCounts[tid] = 0;
    }
    __syncthreads();

//    // copy allele block
//    lAllele[ tid ] = alleles[ allele_idx ];
//    __syncthreads();

    if( tid < BLOCK_PER_STRIDE ) {
        lFree[ tid ] = free_list[ tid ];
    }
    __syncthreads();

    // count free
    unsigned int val = pc.evalGPU( lFree[threadIdx.x] );
    __syncthreads();
    scan_up( val, lCounts );
    // end count free

    if( tid <= BLOCK_PER_STRIDE ) {
        dbg[ tid ] = lCounts[ tid ];
    }
    __syncthreads();
}

__global__ void scatter_mutations( double * rands
                                 , double * alleles
                                 , unsigned int * free_list
                                 , unsigned int * sequences
                                 , unsigned int * mut_events
                                 , unsigned int seq_count
                                 , unsigned int allele_count
                                ) {

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    __shared__ double lMutations[ MAX_THREADS ];                // 8 * 1024
    __shared__ double lAllele[ MAX_THREADS ];                   // 8 * 1024
    __shared__ unsigned int lFree[ BLOCK_PER_STRIDE ];         // 4 * 32
    __shared__ unsigned int lCounts[ BLOCK_PER_STRIDE + 1 ];   // 4 * (32 + 1)
    __shared__ unsigned int lEvents[ MAX_THREADS ];             // 4 * 1024
    __shared__ unsigned int lSeq[ BLOCK_PER_STRIDE ];           // 4 * 32
    __shared__ unsigned int exp_free[ MAX_THREADS ];
                                                                // < 24 K

    popcountGPU< unsigned int > pc;

    unsigned int tmask = (1 << threadIdx.x);

    unsigned int ridx = tid;

    unsigned int sidx = 0, c = 0, d = MAX_THREADS;
    unsigned int allele_idx = tid - MAX_THREADS;    // uses unsigned integer overflow
    unsigned int allele_offset = threadIdx.x - BLOCK_PER_STRIDE;    // uses unsigned integer overflow
    unsigned int block_per_row = allele_count / 32;

    // make sure lCounts is cleared
    if( tid <= BLOCK_PER_STRIDE ) {
        lCounts[tid] = 0;
    }
    __syncthreads();

    while( sidx < seq_count ) {
        if( d == MAX_THREADS ) {
            // load a batch of events into shared memory
            c = (seq_count - sidx + 1);
            if( c > MAX_THREADS ) { c = MAX_THREADS; }

            lEvents[ tid ] = ((tid < c) ? mut_events[ sidx + tid ] : 0);
            c = 0; d = 1;
        }
        __syncthreads();

        unsigned int mStart = lEvents[c], mEnd = lEvents[d];
        c = d++;

        while( mStart < mEnd ) {
            while( lCounts[ BLOCK_PER_STRIDE ] == 0 ) {
                allele_offset += BLOCK_PER_STRIDE;
                allele_idx += MAX_THREADS;

                // copy allele block
                lAllele[ tid ] = alleles[ allele_idx ];
                __syncthreads();

                if( tid < BLOCK_PER_STRIDE ) {
                    lFree[ tid ] = free_list[ allele_offset ];
                }
                __syncthreads();

                // count free
                unsigned int val = pc.evalGPU( lFree[threadIdx.x] );
                __syncthreads();
                scan_up( val, lCounts );
                // end count free
            }
            __syncthreads();

            // number of mutations to be added to sequence
            unsigned int N = (mEnd - mStart);
            if( lCounts[ BLOCK_PER_STRIDE ] < N ) { N = lCounts[ BLOCK_PER_STRIDE ]; }

            if( tid < N ) {
                // only get enough mutations to fill the smaller of either the mutations 
                // for the sequence, or the number of free indices within the stride
                lMutations[ tid ] = rands[ ridx ];
            }
            __syncthreads();

            ridx += N;

            // end count free

            if( tid < BLOCK_PER_STRIDE ) {
                lSeq[ tid ] = sequences[ sidx * block_per_row + allele_offset ];
            }
            __syncthreads();

            unsigned int rstart = lCounts[ threadIdx.y ], rend = lCounts[ threadIdx.y + 1];

            unsigned int free_block = lFree[ threadIdx.y ];
            unsigned int foffset = rstart + pc.evalGPU(free_block & ( tmask - 1 ));
            __syncthreads();

            if( (rstart < rend) && (free_block & tmask) && (foffset < N) ) {
                // replace allele
                lAllele[ tid ] = lMutations[ foffset ];

                // allele no longer free
                exp_free[ tid ] = tmask;
            } else {
                exp_free[ tid ] = 0;
            }
            __syncthreads();

            // reduce exp_free
#pragma unroll
            for( unsigned int i = 2; i <= BLOCK_PER_STRIDE; i <<= 1 ) {
                if( !(tid & (i - 1)) ) {
                    exp_free[tid] |= exp_free[ tid + (i / 2)];
                }
                __syncthreads();
            }

            // not sure if it would be better to use 1 thread from
            // multiple warps to accomplish this task
            //
            // ie. Assuming square thread block <<< WARP, WARP >>>
            // 'simpler' calculation operation but possible
            // uncoalesed memory write as a single thread in each warp
            // performs a write to contiguous memory addcount_ptrses
            //
            // if( threadIdx.x == 0 ) {
            //      lFree[threadIdx.y] &= ~exp_free[ tid ];
            // }
            if( tid < BLOCK_PER_STRIDE ) {
                // update values in shared memory
                unsigned int f = exp_free[ tid * BLOCK_PER_STRIDE ];
                lFree[ tid ] &= ~f;
                lSeq[ tid ] |= f;
            }
            __syncthreads();

            // update free count
            unsigned int v = pc.evalGPU( lFree[ threadIdx.x ] );
            __syncthreads();

            scan_up( v, lCounts );
            // end count free

            // persist modified sequence
            if( tid < BLOCK_PER_STRIDE ) {
                sequences[ sidx * block_per_row + allele_offset ] = lSeq[tid];
            }
            __syncthreads();

            // persist alleles if all free indices have been used
            // next sequence will read alleles 
            if( lCounts[ BLOCK_PER_STRIDE ] == 0 ) {
                alleles[allele_idx] = lAllele[tid];
            }
            __syncthreads();

            mStart += N;
        }

        sidx += 1;
    }

    // copy share back to global
    alleles[ allele_idx ] = lAllele[tid];
    __syncthreads();
}
