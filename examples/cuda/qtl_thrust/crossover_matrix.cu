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

#include "clotho/cuda/warp_sort.hpp"
//#include <iostream>

const unsigned int BLOCK_PER_ROW = 32;
const unsigned int ROW_PER_PAGE = 32;
const unsigned int MAX_THREADS = BLOCK_PER_ROW * ROW_PER_PAGE;

__shared__ double       g_sAlleles[ 1024 ];
__shared__ unsigned int g_sBuffer[ 1024 ];

__constant__ unsigned int g_cEvents[ crossover_wrapper::MAX_EVENTS + 1 ];

//inline std::ostream & operator<<( std::ostream & out, const dim3 & d ) {
//    out << "< " << d.x << ", " << d.y << ", " << d.z << " >";
//    return out;
//}

__global__ void generate_crossover_matrix4( double * rand_pool
                                            , double * alleles
                                            , unsigned int * sequences
                                            , dim3 max_dims ) {
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    unsigned int eStart = g_cEvents[blockIdx.y], eEnd = g_cEvents[blockIdx.y + 1];

    // if there are no recombination events for this thread block/sequence
    if( eStart >= eEnd ) {  // will be true or false for all threads in the block
        if( threadIdx.y == 0 ) {
            sequences[ (blockIdx.y * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.y ] = 0;
        }
        return;
    }

    __shared__ double recomb_points[ 1024 ];
    double rpt = 1.0;

    if( eStart + tid < eEnd ) {
        rpt = rand_pool[ eStart + tid ];
    }
    __syncthreads();

    // order the list of recombination points
    rpt = -log( rpt );

    // scan
    for( unsigned int i = 1; i < 32; i <<= 1 ) {
        double tmp = __shfl_up( rpt, i );
        if( threadIdx.x >= i ) rpt += tmp;
    }

    // share partial sums with other warps in block
    if( threadIdx.x == 31 ) {
        recomb_points[ threadIdx.y ] = rpt;
    }

    double accum = recomb_points[threadIdx.x];
    for( unsigned int i = 1; i < 32; i <<= 1 ) {
        double tmp = __shfl_up( accum, i );
        if( threadIdx.x >= i ) accum += tmp;
    }
    
    if( threadIdx.y > 0 ) {
        double tmp = __shfl( accum, threadIdx.y - 1 );
        rpt += tmp;
    }
    __syncthreads();

    double tmp = __shfl( accum, 31 );
    accum = tmp - log(0.1); // technically should be log of an additional random value (not a constant)

    rpt /= accum;
    recomb_points[tid] = rpt;
    // at this point rpt_{tid} are linearly ordered

    // load an allele into a register
    double allele = alleles[ blockIdx.x * 1024 + tid ];
    __syncthreads();

    // casting to int to allow for negative indices
    int min = 0, max = (eEnd - eStart);

    // binary search of recombination point list
    while( min <= max ) {
        unsigned int mid = ((max - min) / 2) + min;
        
        if( recomb_points[ mid ] < allele ) {
            min = mid + 1;
        } else if( recomb_points[ mid ] > allele ) {
            max = mid - 1;
        } else {
            // allele occurs at a recombination point
            min = mid;
            break;
        }
    }
    __syncthreads();
    // min contains count of preceeding recombination events

    unsigned int mask = ((min & 1) * (1 << threadIdx.x));

    for( unsigned int i = 1; i < 32; i <<= 1) {
        unsigned int m = __shfl_down(mask, i);
        if( !(threadIdx.x & ((i << 1) - 1)) ) mask |= m;
    }

    if(threadIdx.x == 0 ) {
        sequences[ (blockIdx.y * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.y ] = mask;
    }
}

__global__ void generate_crossover_matrix3( double * rand_pool
                                            , double * alleles
                                            , unsigned int * allele_mask
                                            , unsigned int * event_list
                                            , unsigned int * sequences
                                            , dim3 max_dims ) {
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    __shared__ unsigned int sEvents[ 1024 ];
    __shared__ double sOffset[ 1024 ];

    sEvents[ tid ] = event_list[blockIdx.y * 1024 + tid ];
    __syncthreads();

    unsigned int min_events = sEvents[0], max_events = sEvents[ 1023 ];
//    __syncthreads();

    // if there are no recombination events for this sequence
    if( min_events >= max_events ) {  // will be true or false for all threads in the block
        if( threadIdx.y == 0 ) {
            sequences[ blockIdx.y * max_dims.x + (blockIdx.x << 5) + threadIdx.x ] = 0;
        }
        return;
    }

    if( min_events + tid < max_events ) {
        sOffset[ tid ] = rand_pool[ min_events + tid ];
    }
    __syncthreads();    // used to allow for better memory access coalescing in subsequent step

    // load global memory into thread register(s)
    double all = alleles[ blockIdx.x * 1024 + tid ];
    unsigned int all_mask = allele_mask[ blockIdx.x * 1024 + tid ];
     __syncthreads();

    unsigned int idx = all * 1023;  // assumes all in [0,1.0)
    unsigned int lo_event_idx = sEvents[ idx ], hi_event_idx = sEvents[ idx + 1 ];

    lo_event_idx -= min_events;           // number of events before allele bin
    hi_event_idx -= min_events;

    // determine how many events occur **before** current allele in thi_event_idxs bin
    double rng_lo_event_idx = (double) idx / 1024.0, rng_hi_event_idx = ((double)(idx + 1))/ 1024.0;

    // expected to incur branch divergence
    while( lo_event_idx < hi_event_idx ) {
        double val = sOffset[lo_event_idx] * (rng_hi_event_idx - rng_lo_event_idx) + rng_lo_event_idx;
        if( val < all ) {
            // allele in right half of bin
            ++lo_event_idx;
            rng_lo_event_idx = val;
        } else {
            break;
        }
    }

    all_mask *= (lo_event_idx & 1);
    __syncthreads();    // necessary b/c of thread divergence

#pragma unroll
    for( unsigned int i = 2; i <= BLOCK_PER_ROW; i <<= 1 ) {
        unsigned int n = __shfl_down(all_mask, (i >> 1), BLOCK_PER_ROW );
        if( !(threadIdx.x & (i - 1))) all_mask |= n;
    }

    if( threadIdx.x == 0 ) {
        sequences[ blockIdx.y * max_dims.x + (blockIdx.x << 5) + threadIdx.y ] = all_mask;
    }
}

// kernel is invoked with each sequence as its own block Y coordinate
__global__ void generate_crossover_matrix2a( double * rand_pool
                                            , double * alleles
                                            , unsigned int * event_list
                                            , unsigned int * sequences
                                            , dim3 max_dims ) {
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    __shared__ double sRec[ 1024 ];

    unsigned int eStart = event_list[ max_dims.y +  blockIdx.y ], eEnd = event_list[ max_dims.y +  blockIdx.y + 1 ];
    __syncthreads();

    if( eStart >= eEnd ) {  // should be true or false for all threads in the block
        if( threadIdx.y == 0 ) {
            sequences[ blockIdx.y * max_dims.x + (blockIdx.x << 5) + threadIdx.x ] = 0;
        }
        __syncthreads();
        return;
    }

    double all = alleles[ blockIdx.x * 1024 + tid ];
    __syncthreads();

    unsigned int N = eEnd - eStart;
    if( tid < N ) {
        sRec[tid] = rand_pool[ eStart + tid ];
    }
    __syncthreads();

    double rng_lo = 0.0, rng_hi = 1.0;
    unsigned int prior_events = 0, event_idx = 0;
    while( N-- ) {
        double offset = sRec[ event_idx ];
        double r = rng_lo + offset * (rng_hi - rng_lo);

        unsigned int n = offset * N;
        if( r < all ) {
            // allele in right half of bin
            prior_events += n + 1;
            event_idx += n + 1;
            N -= n;
            rng_lo = r;
        } else {
            ++event_idx;
            N = n;
            rng_hi = r;
        }
    }
    __syncthreads();

    // if there are an odd number of preceeding crossover events
    // then this allele should be selected from the 'other' sequence
    // (or it exists in a crossover region)
    unsigned int m = ((prior_events & 1) * (1 << threadIdx.x));

#pragma unroll
    for( unsigned int i = 2; i <= BLOCK_PER_ROW; i <<= 1 ) {
        unsigned int n = __shfl_down(m, (i >> 1), BLOCK_PER_ROW );
        if( !(threadIdx.x & (i - 1))) m |= n;
    }

    // use one thread in each warp to write back to global memory
    // expecting this to not coalesce well
    // 
    if( threadIdx.x == 0 ) {
        sequences[ blockIdx.y * max_dims.x + (blockIdx.x << 5) + threadIdx.y ] = m;
    }
}

// kernel is invoked with each sequence as its own block Y coordinate
__global__ void generate_crossover_matrix2( double * rand_pool
                                            , double * alleles
                                            , unsigned int * event_list
                                            , unsigned int * sequences
                                            , dim3 max_dims ) {
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
//    unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;

//    __shared__ unsigned int sMask[ 32 ];
    __shared__ double sRec[ 1024 ];

    unsigned int eStart = event_list[ max_dims.y +  blockIdx.y ], eEnd = event_list[ max_dims.y +  blockIdx.y + 1 ];
    __syncthreads();

    if( eStart >= eEnd ) {  // should be true or false for all threads in the block
        if( threadIdx.y == 0 ) {
            sequences[ blockIdx.y * max_dims.x + (blockIdx.x << 5) + threadIdx.x ] = 0;
        }
        __syncthreads();
        return;
    }

    double all = alleles[ blockIdx.x * 1024 + tid ];
    __syncthreads();

    unsigned int mask = (1 << threadIdx.x);
    unsigned int m = 0;

    while( eStart < eEnd ) {
        __syncthreads();

        unsigned int N = (eEnd - eStart);
        if( N > 1024 ) { N = 1024; }

        sRec[ tid ] = (( tid < N ) ? rand_pool[ eStart + tid ] : -1.0);
        __syncthreads();

        eStart += N;

        while( N-- ) {
            double r = sRec[ N ];
            __syncthreads();
            m ^= (( 0 <= all && all < r) ? mask : 0 );
        }
    }
    __syncthreads();

#pragma unroll
    for( unsigned int i = 2; i <= BLOCK_PER_ROW; i <<= 1 ) {
        unsigned int n = __shfl_down(m, (i >> 1), BLOCK_PER_ROW );
        if( !(threadIdx.x & (i - 1))) m |= n;
    }
//    if( threadIdx.x == 0 ) {
//        sMask[ threadIdx.y ] = m;
//    }
//    __syncthreads();
//
//    if( threadIdx.y == 0 ) {
//        sequences[ blockIdx.y * max_dims.x + blockIdx.x * 32 + threadIdx.x ] = sMask[threadIdx.x];
//    }

    // use one thread in each warp to write back to global memory
    // expecting this to not coalesce well
    // 
    if( threadIdx.x == 0 ) {
        sequences[ blockIdx.y * max_dims.x + (blockIdx.x << 5) + threadIdx.y ] = m;
    }
}

__global__ void generate_crossover_matrix( double * rand_pool
                                         , double * alleles
                                         , unsigned int * event_list
                                         , unsigned int * sequences
                                         , dim3 max_dims ) {
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    // all threads within a block will operate upon a contiguous set of events and alleles
    __shared__ double lEvents[ MAX_THREADS ];       //  8K
    __shared__ unsigned int lMask[ 32 ];            //  4K
    __shared__ unsigned int lCounts[ MAX_THREADS ]; //  4K
                                                    // 24K shared memory

    // copy the alleles that this block is going to be inspecting
    unsigned int aidx = blockIdx.x * MAX_THREADS + tid;
    double all = ((aidx < max_dims.x ) ? alleles[aidx] : -1.0);
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

        // number of events for sequence
        unsigned int eStart = lCounts[ c ], eEnd = lCounts[ d ];
        __syncthreads();
        c = d++;

        if( eStart < eEnd ) {
            unsigned int m = 0;
            do {
                // refresh event list

                // current method performs a 'lookahead copy'.
                // general question: is it better to perform a block copy from global memory to local memory,
                // or to simply access global memory as necessary?

                // Basic concern is that when N is small there will be a lot of sleeping threads
                unsigned int N = (eEnd - eStart);
                if( N > MAX_THREADS ) { N = MAX_THREADS; }
                lEvents[ tid ] = ((tid < N ) ? rand_pool[ eStart + tid ] : 0);
            
                __syncthreads();

                eStart += N;

                while( N-- ) {
                    double ev = lEvents[ N ];
                    __syncthreads();

                    m ^= (( 0 <= all && all < ev) ? mask : 0 );
                }
            } while( eStart < eEnd );
           __syncthreads();

#pragma unroll
            for( unsigned int i = 2; i <= BLOCK_PER_ROW; i <<= 1 ) {
                unsigned int n = __shfl_down(m, (i / 2), BLOCK_PER_ROW );
                if( !(threadIdx.x & (i - 1))) m |= n;
            }

            if( threadIdx.x == 0 ) {
                lMask[ threadIdx.y ] = m;
            }
        } else {
            if( threadIdx.y == 0 ) {
                lMask[ threadIdx.x ] = 0;
            }
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
    g_sAlleles[ tid ] = ((tid < count ) ? alleles[ tid ] : -1.0);
}

__global__ void init_sequence( unsigned int * seq ) {
    unsigned int tid = threadIdx.y* blockDim.x + threadIdx.x;
    if( threadIdx.y == 0 ) {
        g_sBuffer[ threadIdx.x ] = seq[ threadIdx.x ];
    }
    __syncthreads();

    unsigned int mask = (1 << threadIdx.x );
    unsigned int v = g_sBuffer[ threadIdx.x ];
    __syncthreads();

    g_sBuffer[ tid ] = (v & mask);
}

__global__ void crossover2( double rpoint) {

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    double a = g_sAlleles[ tid ];
    unsigned int res = g_sBuffer[ tid ];
    __syncthreads();

    unsigned int mask = (1 << threadIdx.x);
    unsigned int b = !!((0.0 <= a) && (a < rpoint)); // rec point is after allele position
    res = (res ^ ( b*mask));

    g_sBuffer[tid] = res;
}

__global__ void finalize_sequence( unsigned int * seq ) {
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    unsigned int v = g_sBuffer[tid];
    __syncthreads();

#pragma unroll
    for( unsigned int i = 2; i <= 32; i<<=1 ) {
        unsigned int tmp = __shfl_down( v, (i / 2), 32 );
        if( !(threadIdx.x & (i - 1)) ) v |= tmp;
    }

    if( threadIdx.x == 0 ) {
        g_sBuffer[ threadIdx.y ] = v;
    }
    __syncthreads();

    if( threadIdx.y == 0 ) {
        seq[threadIdx.x] = g_sBuffer[threadIdx.x];
    }
}

void crossover_wrapper::operator()( double * rand_pool, double * alleles, unsigned int * event_list, unsigned int * sequences, dim3 max_dims ) {

    unsigned int seq_count = max_dims.y;
    while( seq_count ) {
        unsigned int N = seq_count;
        if( N > MAX_EVENTS ) { N = MAX_EVENTS; }
        cudaMemcpyToSymbol( event_list, g_cEvents, N + 1, cudaMemcpyDeviceToDevice );

        dim3 blocks( max_dims.x / 1024, N, 1), threads( 32, 32, 1);

        generate_crossover_matrix4<<< blocks, threads >>>( rand_pool, alleles, sequences, max_dims );


//        std::cerr << blocks << ", " << threads << std::endl;
        cudaDeviceSynchronize();

        event_list += N;
        seq_count -= N;
        sequences += (blocks.y * blocks.x * threads.x);
    }
}
