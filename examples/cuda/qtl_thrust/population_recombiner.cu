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
#include "population_recombiner.hpp"

const unsigned int SEQ_PER_INDIVIDUAL = 2;
const unsigned int PARENT_PER_OFFSPRING = 2;
const unsigned int RANDOM_PER_PARENT = 2;   // id offset & swapping of source sequences

const unsigned int MAX_OFFSPRING = 8;
const unsigned int MAX_BLOCKS_PER_STRIDE = 32;

#define compute_randoms_per_offspring(x) RANDOM_PER_PARENT * PARENT_PER_OFFSPRING * x

const unsigned int MAX_OFFSPRING_SEQ = SEQ_PER_INDIVIDUAL * MAX_OFFSPRING;  // 16
const unsigned int MAX_PARENTS = PARENT_PER_OFFSPRING * MAX_OFFSPRING;      // 16
const unsigned int MAX_PARENT_SEQ = SEQ_PER_INDIVIDUAL * MAX_PARENTS;       // 32

const unsigned int MAX_PARENT_BLOCKS = MAX_PARENT_SEQ * MAX_BLOCKS_PER_STRIDE;  // 32 * 32 = 1024
const unsigned int MAX_OFFSPRING_BLOCKS = MAX_OFFSPRING_SEQ * MAX_BLOCKS_PER_STRIDE;    // 16 * 32 = 512
const unsigned int MAX_RANDOMS = compute_randoms_per_offspring( MAX_OFFSPRING );   // 2 * 2 * 8 = 32

const unsigned int MAX_THREADS = MAX_OFFSPRING_BLOCKS;  // 512

__device__ unsigned int get_parent_sequence_index( double offset, unsigned int tot_seqs) {
    // uses integer division to obtain the floor sequence index
    unsigned int id = offset * (tot_seqs / SEQ_PER_INDIVIDUAL);
    id = SEQ_PER_INDIVIDUAL * id;
    return id;
}

/*
 * Is it better to have 1024 threads to perform a single round of memory loading,
 * only to turn around an use half (or less) to perform computation? Or, should
 * only 512 threads be allocated, perform two loads, then a single reduction?
 *
 * Opted for the double load approach because it simplifies logic
 */
__global__ void recombine_population( double * rand_pool
                                    , unsigned int * parents
                                    , unsigned int * offspring
                                    , unsigned int parent_rows
                                    , unsigned int parent_cols
                                    , unsigned int off_rows
                                    , unsigned int off_cols) {
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    __shared__ double lRandom[ MAX_RANDOMS ];
    __shared__ unsigned int lParents[ MAX_PARENT_BLOCKS ];
    __shared__ unsigned int lOffspring[ MAX_OFFSPRING_BLOCKS ];

    unsigned int off_per_block = off_rows / blockDim.y;

    unsigned int off_start = threadIdx.y * off_per_block, off_end = off_start + off_per_block;

    if( off_end > off_cols ) { off_end = off_cols; }

    unsigned int boff = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int rand_idx =  RANDOM_PER_PARENT * off_start + tid;

    unsigned int N = (off_end - off_start) / MAX_OFFSPRING;

    while( N-- ) {
        if( tid < MAX_RANDOMS ) {
            lRandom[ tid ] = rand_pool[ rand_idx ];
            rand_idx += MAX_RANDOMS;
        }

        __syncthreads();

        unsigned int idx = get_parent_sequence_index( lRandom[ (threadIdx.y & 0xFFFFFFFE) ], parent_rows);
        idx = idx * parent_cols + boff;
        // copy parents
        lParents[ tid ] = parents[ idx ];   // parent sequence

        idx += parent_cols;
        lParents[ tid + MAX_THREADS ] = parents[ idx ];   // parent sequence

        // copy crossovers
        idx = off_start * off_cols + boff;
        lOffspring[ tid ] = offspring[ idx ];
        __syncthreads();

        unsigned int cross_mask = lOffspring[ tid ];
        if( lRandom[ ((threadIdx.y & 0xFFFFFFFE) + 1) ] < 0.5 ) {
            cross_mask = ~cross_mask;
        }
        lParents[ tid ] = ((lParents[ tid ] & cross_mask) | (lParents[ tid + MAX_THREADS ] & ~cross_mask));
        __syncthreads();

        // copy sequence back to offspring
        offspring[ idx ] = lParents[ tid ];
        __syncthreads();

        off_start += MAX_OFFSPRING;
    }

    __syncthreads();

    // handling tail separately to simplify logic
    // from coding perspective, this produces duplicate code which is bad
    // should functionalize the inner logic
    if( off_start < off_end ) {
        N = off_end - off_start;    // offspring left to build
        unsigned int seqs = compute_randoms_per_offspring( N ); // 2 * 2 * N
        if( tid < seqs ) {
            lRandom[ tid ] = rand_pool[ rand_idx ];
        }

        __syncthreads();

        unsigned int offidx = off_start * off_cols + boff;

        // since the blockDim.x is assumed to be the size of warp
        // these if blocks will allow specific warps to be skip execution
        // rather than specific threads within a warp
        bool use_thread = (threadIdx.y < 2 * N);
        if( use_thread ) { // for each parent
            unsigned int idx = get_parent_sequence_index( lRandom[ (threadIdx.y & 0xFFFFFFFE) ], parent_rows);
            idx = idx * parent_cols + boff;

            lParents[ tid ] = parents[ idx ];   // parent first sequence

            idx += parent_cols;
            lParents[ tid + MAX_THREADS ] = parents[ idx ];   // parent second sequence

            lOffspring[ tid ] = offspring[ offidx ];
        } else {
            lParents[ tid ] = 0;    // clear the rest
            lParents[ tid + MAX_THREADS ] = 0;
            lOffspring[tid] = 0;
        }
        __syncthreads();

        if( use_thread ) {
            unsigned int cross_mask = lOffspring[ tid ];
            if( lRandom[ ((threadIdx.y & 0xFFFFFFFE) + 1) ] < 0.5 ) {
                cross_mask = ~cross_mask;
            }
            lParents[ tid ] = ((lParents[ tid ] & cross_mask) | (lParents[ tid + MAX_THREADS ] & ~cross_mask));
        }
        __syncthreads();

        // copy sequence back to offspring
        if( use_thread ) {
            offspring[ offidx ] = lParents[ tid ];
        }
        __syncthreads();
    }
}

__global__ void recombine( unsigned int * p0
                            , unsigned int * p1
                            , unsigned int * off
                            , unsigned int cols ) {
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    unsigned int boffset = blockIdx.x * blockDim.x + tid;

    unsigned int p = ((boffset < cols) ? p0[ boffset ] : 0 );
    unsigned int q = ((boffset < cols) ? p1[ boffset ] : 0 );
    unsigned int res = ((boffset < cols) ? off[ boffset ] : 0 );
    __syncthreads();

    res = (( p & ~res ) | ( q & res ));
    __syncthreads();

    if( boffset < cols ) {
        off[ boffset ] = res;
    }
}

__global__ void recombiner( double * rands
                            , unsigned int * parents
                            , unsigned int parent_rows
                            , unsigned int parent_cols
                            , unsigned int * off
                            , unsigned int cols
                            , unsigned int seq_offset ) {
    double id_offset = rands[ seq_offset + blockIdx.y ];
    __syncthreads();

    unsigned int col_offset = (blockIdx.x + threadIdx.y) * blockDim.x + threadIdx.x;

    // using integer cast to truncate of fractional portion
    unsigned int p0_offset = id_offset * ((parent_rows - 1) / 2);
    p0_offset = (2 * p0_offset * parent_cols) + col_offset;

    unsigned int p = 0, q = 0, res = 0;
    if( col_offset < parent_cols ) {
        // should hold true for entire warps
        p = parents[ p0_offset ];
        q = parents[ p0_offset + parent_cols ];
    }
    __syncthreads();

    if( col_offset < cols ) {
        res = off[ (seq_offset + blockIdx.y) * cols + col_offset ];
    }
    __syncthreads();

    res = (( p & ~res ) | ( q & res ));
    __syncthreads();

    if( col_offset < cols ) {
        off[ (seq_offset + blockIdx.y) * cols + col_offset ] = res;
    }
}
