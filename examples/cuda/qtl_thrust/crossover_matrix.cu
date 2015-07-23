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

const unsigned int BLOCK_PER_ROW = 32;
const unsigned int ROW_PER_PAGE = 32;
const unsigned int MAX_THREADS = BLOCK_PER_ROW * ROW_PER_PAGE;

__shared__ double f_sAlleles[ 1024 ];

__global__ void generate_crossover_matrix( double * rand_pool
                                         , double * alleles
                                         , unsigned int * event_list
                                         , unsigned int * sequences
                                         , dim3 max_dims ) {
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    // all threads within a block will operate upon a contiguous set of events and alleles
    __shared__ double lEvents[ MAX_THREADS ];       //  8K
    __shared__ double lAlleles[ MAX_THREADS ];      //  8K
    __shared__ unsigned int lMask[ MAX_THREADS ];   //  4K
    __shared__ unsigned int lCounts[ MAX_THREADS ]; //  4K
                                                    // 24K shared memory

    // copy the alleles that this block is going to be inspecting
    unsigned int aidx = blockIdx.x * MAX_THREADS + tid;
    lAlleles[ tid ] = ((aidx < max_dims.x ) ? alleles[aidx] : -1.0);
    __syncthreads();

    unsigned int seq_per_block = max_dims.y / gridDim.y;
    unsigned int unit_per_seq = max_dims.x / 32;

    unsigned int seq_start = blockIdx.y * seq_per_block, seq_end = seq_start + seq_per_block;

    if( seq_end > max_dims.y ) seq_end = max_dims.y;

    unsigned int c = 0, d = MAX_THREADS;
    unsigned int mask = (1 << threadIdx.x);

    while( seq_start < seq_end ) {
        if( d == MAX_THREADS ) {
            // refresh event count list
            c = seq_end - seq_start + 1;
            if( c > MAX_THREADS ) { c = MAX_THREADS; }

            lCounts[ tid ] = (( tid < c ) ? event_list[ seq_start + tid ] : 0);
            c = 0;
            d = 1;
        }
        __syncthreads();

        lMask[tid] = 0; // reset masks for sequence
        __syncthreads();

        // number of events for sequence
        unsigned int eStart = lCounts[ c ], eEnd = lCounts[ d ];
        c = d++;

        while( eStart < eEnd ) {
            // refresh event list

            // current method performs a 'lookahead copy'.
            // general question: is it better to perform a block copy from global memory to local memory,
            // or to simply access global memory as necessary?

            // Basic concern is that when N is small there will be a lot of sleeping threads
            unsigned int N = eEnd - eStart;
            if( N > MAX_THREADS ) { N = MAX_THREADS; }
            lEvents[ tid ] = ((tid < N ) ? rand_pool[ eStart + tid ] : 0);
            
            __syncthreads();

            eStart += N;

            while( N-- ) {
                lMask[tid] ^= (( lAlleles[tid] >= 0 && lAlleles[tid] < lEvents[N]) ? mask : 0 );
            }
        }

        __syncthreads();
/*
        // collapse masks for sequence
        if( !(tid & 1) ) {
            // half warp
            // merge odds and evens: 0 < 1; 2 < 3; 4 < 5 ...
            lMask[tid] |= lMask[tid + 1];
        }

        __syncthreads();

        if( !(tid & 3) ) {
            // quarter warp
            // merge 0 < 2;4 < 6;8 < 10...
            lMask[tid] |= lMask[tid + 2];
        }

        __syncthreads();

        if( !(tid & 7) ) {
            // eighth warp
            // merge 0 < 4; 8 < 12; 16 < 20; 24 < 28
            lMask[tid] |= lMask[tid + 4];
        }

        __syncthreads();

        if(! (tid & 15) ) {
            // sixth warp
            // merge 0 < 8; 16 < 24
            lMask[tid] |= lMask[tid + 8];
        }

        __syncthreads();

        // copy sequence to output
        if( !(tid & 31) ) {
            // thirty-second warp (1-thread)
            // merge 0 < 16
            lMask[ tid ] |=  lMask[ tid + 16 ];
        }

        __syncthreads();*/

        unsigned int m = lMask[tid];
#pragma unroll
        for( unsigned int i = 2; i <= BLOCK_PER_ROW; i <<= 1 ) {
            unsigned int n = __shfl_down(m, (i / 2), BLOCK_PER_ROW );
            if( !(threadIdx.x & (i - 1))) m |= n;
        }

        if( threadIdx.x == 0 ) {
            lMask[ threadIdx.y ] = m;
        }
        __syncthreads();

        if( threadIdx.y == 0 ) { // threads [0, BLOCK_PER_ROW )
            // use a single warp to copy between shared and global memory
            sequences[ seq_start * unit_per_seq + tid ] =  lMask[threadIdx.x];
        }
        
        __syncthreads();

        ++seq_start;
    }
}

/*
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

    unsigned int seq_start = blockIdx.y * seq_per_block, seq_end = seq_start + seq_per_block;

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
            sequences[ seq_start * max_dims.x / 32 + tid ] = lMask[ tid * 32 ];
        }
        
        __syncthreads();
        ++seq_start;
    }
}*/

__global__ void crossover( unsigned int * seq
                            , double * alleles
                            , unsigned int allele_count
                            , double rpoint ) {
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    __shared__ unsigned int sSeq[32];

    // Thought: every thread block copies a relative
    // tile of the global allele and free_list
    // Technically these could exist in
    // shared memory as another kernel call with
    // the same blockIdx coordinates will copy
    // the same relative tile back into shared memory
    // ! Need to determine whether this repeated global access
    // is a significant overhead.
    // !! From an algorithmic perspective this may
    // be most general approach as it would avoid
    // scenarios where shared memory is exhausted
    unsigned int aidx = 32 * blockIdx.x * blockDim.x + tid;
    double sAllele = ((aidx < allele_count) ? alleles[aidx] : -1.0);
    __syncthreads();

    // load registers with global value
    unsigned int s = seq[ blockIdx.x * blockDim.x + threadIdx.x ];
    __syncthreads();

    unsigned int mask = (1 << threadIdx.x );
    unsigned int b = !!( (sAllele >= 0.0) && (sAllele < rpoint)); // rec point is after allele position
    s = ((s & mask) ^ ( b*mask));

    __syncthreads();

#pragma unroll
    for( unsigned int i = 2; i <= 32; i<<=1 ) {
        unsigned int tmp = __shfl_down( s, (i / 2), 32 );
        if( !(threadIdx.x & (i - 1)) ) s |= tmp;
    }

    if( threadIdx.x == 0 ) {
        sSeq[ threadIdx.y ] = s;
    }
    __syncthreads();

    if( threadIdx.y == 0 ) {
        seq[ blockIdx.x * blockDim.x + threadIdx.y ] = sSeq[threadIdx.x];
    }
}

__global__ void init_alleles( double * alleles, unsigned int count ) {
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    f_sAlleles[ tid ] = ((tid < count ) ? alleles[ tid ] : -1.0);
}

__global__ void crossover2( unsigned int * seq, double rpoint) {

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    double a = f_sAlleles[ tid ];
    __syncthreads();

    __shared__ unsigned int sSeq[ 32 ];
    if( threadIdx.y == 0 ) {
        sSeq[ threadIdx.x ] = seq[ threadIdx.x ];
    }
    __syncthreads();

    unsigned int res = sSeq[ threadIdx.x ];
    __syncthreads();

    
    unsigned int mask = (1 << threadIdx.x );
    unsigned int b = !!((0.0 <= a) && (a < rpoint)); // rec point is after allele position
    res = ((res & mask) ^ ( b*mask));

#pragma unroll
    for( unsigned int i = 2; i <= 32; i<<=1 ) {
        unsigned int tmp = __shfl_down( res, (i / 2), 32 );
        if( !(threadIdx.x & (i - 1)) ) res |= tmp;
    }

    if( threadIdx.x == 0 ) {
        sSeq[ threadIdx.y ] = res;
    }
    __syncthreads();

    if( threadIdx.y == 0 ) {
        seq[threadIdx.x] = sSeq[threadIdx.x];
    }
}
