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
#ifndef CROSSOVER_KERNEL_UNORDERED_IMPL_HPP_
#define CROSSOVER_KERNEL_UNORDERED_IMPL_HPP_

#include "clotho/cuda/crossover/crossover_config_def.hpp"

#include "clotho/cuda/data_spaces/allele_space/device_allele_space.hpp"
#include "clotho/cuda/data_spaces/sequence_space/device_sequence_space.hpp"
#include "clotho/cuda/data_spaces/free_space/device_free_space.hpp"
#include "clotho/cuda/distributions/poisson_distribution.hpp"

#include "clotho/cuda/data_spaces/tags/unordered_tag.hpp"
#include "clotho/cuda/crossover/persist_sequence.hpp"
#include "clotho/utility/algorithm_version.hpp"


/**
 * Assumptions:
 *  - 1 sequence/tile
 *  - 1 tile/block
 *  - 4 warp/tile (1 thread/allele * 128 allele/ALIGNMENT_SIZE) * 1/32 warp/thread * 1 ALIGNMENT_SIZE/tile
 */
template < class StateType, class AlleleSpaceType, class RealType, class IntType, unsigned char V >
__global__ void crossover_kernel( StateType * states
                                , AlleleSpaceType * alleles
                                , device_free_space< IntType, unordered_tag > * free_space
                                , poisson_cdf< RealType, 32 > * pois
                                , device_sequence_space< IntType > * sequences
                                , clotho::utility::algo_version< V > * v ) {

    typedef StateType                                   state_type;
    typedef AlleleSpaceType                             allele_space_type;
    typedef typename allele_space_type::real_type       real_type;

    typedef device_sequence_space< IntType >            sequence_space_type;
    typedef typename sequence_space_type::int_type      int_type;

    typedef poisson_cdf< RealType, 32 >                 poisson_type;
    typedef xover_config< unordered_tag, V >            xover_type;

    assert( blockDim.x == 32 && blockDim.y <= xover_type::MAX_WARPS );

    int_type bid = blockIdx.y * gridDim.x + blockIdx.x;
    int_type tid = threadIdx.y * blockDim.x + threadIdx.x;

    unsigned int tpb = (blockDim.x * blockDim.y);
    assert( (tpb % allele_space_type::ALIGNMENT_SIZE) == 0 );

    unsigned int bpg = (gridDim.x * gridDim.y);

    unsigned int i;
    int_type sequence_width = sequences->seq_width;
    int_type nSequences = sequences->seq_count;    

    int_type    * seqs = sequences->sequences;

    real_type       * allele_list = alleles->locations;
    unsigned int    nAlleles = alleles->capacity;

    assert( (nAlleles % allele_space_type::ALIGNMENT_SIZE) == 0 );

    if( nAlleles == 0 ) { return; }

    __shared__ real_type        s_pois_cdf[ poisson_type::MAX_K ];
    __shared__ real_type        rand_pool[ allele_space_type::ALIGNMENT_SIZE ];
    __shared__ unsigned int     event_hash[ allele_space_type::ALIGNMENT_SIZE];

    unsigned int max_k = pois->max_k;
    if( tid < poisson_type::MAX_K ) {
        s_pois_cdf[ tid ] = pois->_cdf[tid];
    }
    __syncthreads();

    state_type local_state = states[ bid * tpb + tid ];

    unsigned int nonzero_warp = (threadIdx.y != 0);
    unsigned int nonzero_thread = (tid != 0);
    int_type seq_idx = bid;
    while( seq_idx < nSequences ) {
        real_type x = curand_uniform( &local_state );
        rand_pool[ tid ] = curand_uniform( &local_state );

        int_type rand = _find_poisson_maxk32( s_pois_cdf, x, max_k );
        __syncthreads();

        for( i = 1; i < 32; i <<= 1 ) {
            unsigned int r = __shfl_up( rand, i );
            rand += ( (threadIdx.x >= i ) * r );
        }

        if( threadIdx.x == 31 ) {
            event_hash[ threadIdx.y ] = rand;
        }
        __syncthreads();

        unsigned int _sum = event_hash[ threadIdx.x ];
        _sum *= (threadIdx.x < blockDim.y);
        __syncthreads();

        for( i = 1; i < 32; i <<= 1 ) {
            unsigned int s = __shfl_up( _sum, i );
            _sum += (( threadIdx.x >= i ) * s);
        }

        unsigned int s = __shfl( _sum, 31 );
//        assert( max_events < allele_space_type::ALIGNMENT_SIZE );

        s = __shfl( _sum, threadIdx.y - nonzero_warp);
        s *= nonzero_warp;
        __syncthreads();

        rand += s;
        event_hash[tid] = rand;
        __syncthreads();

        i = event_hash[ tid - nonzero_thread];    // minimum event index
        i *= nonzero_thread;
        __syncthreads();

        // BEGIN divergent code
        real_type accum = 0.;
        while (i < rand) {
            x = rand_pool[ i ];

            accum += (log( x ) / (real_type)(rand - i));

            rand_pool[i++] = ((((real_type)tid) + (1.0 - exp(accum))) / ((real_type)allele_space_type::ALIGNMENT_SIZE));
        }
        __syncthreads();
        // END divergent code

        unsigned int seq_offset = seq_idx * sequence_width + threadIdx.y;   // every thread in a warp has the same block offset
        i = tid;
        while( i < nAlleles ) {
            x = allele_list[ i ];

            rand = (unsigned int) ( x * ((real_type)allele_space_type::ALIGNMENT_SIZE));

            unsigned int nonzero_bin = (rand != 0);  // _keep == 0 -> rand == 0 -> e_min == 0; _keep == 1 -> rand != 0 -> e_min == event_hash[ rand - 1]

            // each thread reads hash (and random pool) relative to their 
            // local allele (x)
            // this code will result in bank conflicts
            //
            // initial performance results suggest that this is
            // an acceptable overhead as overall runtime of simulation
            // loop is minimized when this algorithm is used
            unsigned int e_max = event_hash[ rand ];
            unsigned int e_min = event_hash[ rand - nonzero_bin ];
            e_min *= nonzero_bin;

            int_type cmask = e_min;

            // BEGIN divergent code
            while( e_min < e_max ) {
                accum = rand_pool[ e_min++ ];
                cmask += (accum < x);
            }
            __syncthreads();
            // END divergent code

            cmask = ((cmask & 1) << threadIdx.x);
            
            for( unsigned int j = 1; j < 32; j <<= 1 ) {
                int_type _c = __shfl_up( cmask, j);
                cmask |= ((threadIdx.x >= j) * _c);
            }
            if( threadIdx.x == 31 ) {
                seqs[ seq_offset ] = cmask;
            }
            __syncthreads();

            i += tpb;
            seq_offset += blockDim.y;   // wpb
        }
        __syncthreads();


        seq_idx += bpg;
    }

    states[ bid * tpb + tid ] = local_state;
}

template < class StateType, class AlleleSpaceType, class RealType, class IntType >
__global__ void crossover_kernel( StateType * states
                                , AlleleSpaceType * alleles
                                , device_free_space< IntType, unordered_tag > * free_space
                                , poisson_cdf< RealType, 32 > * pois
                                , device_sequence_space< IntType > * sequences
                                , clotho::utility::algo_version< 2 > * v ) {

    typedef StateType  state_type;
    typedef AlleleSpaceType allele_space_type;
    typedef typename allele_space_type::real_type   real_type;

    typedef device_sequence_space< IntType >            sequence_space_type;
    typedef typename sequence_space_type::int_type      int_type;

    typedef poisson_cdf< RealType, 32 > poisson_type;

    unsigned int  nAlleles = alleles->capacity;
    if( nAlleles == 0 ) { return; }

    assert( (nAlleles & 31) == 0 ); // multiple of 32 alleles

    const unsigned int MAX_EVENTS_PER_WARP = 64;    // maximum number of recombination events per sequence
    const unsigned int MAX_WARP_PER_BLOCK = 32;
    const unsigned int HASH_WIDTH = 32;

    unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int lane_id = (tid & 31);
    unsigned int warp_id = (tid >> 5);

    unsigned int tpb = (blockDim.x * blockDim.y);
    assert( (tpb & 31) == 0 );  // all warps are full

    unsigned int bpg = (gridDim.x * gridDim.y);
    unsigned int wpb = (tpb >> 5);  // == sequences/block (spb)

    assert( wpb <= MAX_WARP_PER_BLOCK );

    unsigned int spg = bpg * wpb;
    
    int_type sequence_width = sequences->seq_width;
    int_type nSequences = sequences->seq_count;    

    int_type    * seqs = sequences->sequences;
    real_type   * allele_list = alleles->locations;

    __shared__ real_type        s_pois_cdf[ poisson_type::MAX_K ];                      // 4 * 32 == 128
    __shared__ real_type        rand_pool[ MAX_WARP_PER_BLOCK * MAX_EVENTS_PER_WARP ];  // 4 * 32 * 64 == 8192
    __shared__ unsigned int     event_hash[ MAX_WARP_PER_BLOCK * HASH_WIDTH ];          // 4 * 32 * 128 == 16384
                                                                                        //----------------------
                                                                                        // 24704 (bytes of shared memory)

    unsigned int max_k = pois->max_k;
    if( tid < poisson_type::MAX_K ) {
        s_pois_cdf[ tid ] = pois->_cdf[tid];
    }
    __syncthreads();

    state_type local_state = states[ bid * tpb + tid ]; // assumes every thread in the GRID has a state defined

    unsigned int max_seq_id = nSequences / wpb;
    max_seq_id += (( nSequences % wpb ) ? 1 : 0);
    max_seq_id *= wpb;

    int_type seq_idx = bid * wpb + warp_id;

    while( seq_idx < max_seq_id ) { // should allow blocks to terminate early

        // generate a recombination hash for each sequence (warp)
//        unsigned int psum = 0;
//        for( unsigned int w = lane_id; w < HASH_WIDTH; w += 32 ) {
            real_type x = curand_uniform( &local_state );

            unsigned int rand = _find_poisson_maxk32( s_pois_cdf, x, max_k );
            __syncthreads();

            // compute prefix sum with in each warp
            for( unsigned int i = 1; i < 32; i <<= 1 ) {
                unsigned int r = __shfl_up( rand, i );
                rand += ( (lane_id >= i ) * r );
            }

//            rand += psum;
//            event_hash[ warp_id * HASH_WIDTH + w ] = rand;   //  event_hash contain prefix sum (scan)
            event_hash[ tid ] = rand;

            unsigned int s = __shfl_up( rand, 1);   //
            s *= ( lane_id != 0 );

            s += (warp_id * MAX_EVENTS_PER_WARP);       // shift s and rand to be relative to warp
            rand += (warp_id * MAX_EVENTS_PER_WARP);

            // BEGIN divergent code
            real_type accum = 0.;
            while (s < rand) {
                x = curand_uniform( &local_state );

                accum += (log( x ) / (real_type)(rand - s));

                rand_pool[s++] = ((((real_type)tid) + (1.0 - exp(accum))) / ((real_type)32.));
            }
            __syncthreads();
            // END divergent code
            //
//            psum = __shfl( rand, 31 );
//        }

        unsigned int seq_offset = seq_idx * sequence_width;
        unsigned int a_id = lane_id;
        while( a_id < nAlleles ) {
            real_type x = allele_list[ a_id ];

            unsigned int h_idx = (unsigned int) ( x * ((real_type) HASH_WIDTH ));  // map allele location to bin index

            unsigned int s = event_hash[ h_idx++ ];
            unsigned int e = event_hash[ h_idx ];
            __syncthreads();

            int_type cmask = s;

            s += (warp_id * MAX_EVENTS_PER_WARP);
            e += (warp_id * MAX_EVENTS_PER_WARP );

            // BEGIN divergent code
            while( s < e  ) {
                real_type y = rand_pool[ s++ ];
                cmask += (x > y);
            }
            __syncthreads();
            // END divergent code

            cmask = ((cmask & 1) << lane_id);

            // reduce cmask within each warp
            for( unsigned int i = 1; i < 32; i <<= 1 ) {
                int_type _c = __shfl_up( cmask, i );
                cmask |= ((lane_id >= i) * _c);
            }

            if( lane_id == 31 && seq_idx < nSequences ) {
                seqs[seq_offset] = cmask;
            }
            __syncthreads();

            a_id += 32;
            ++seq_offset;
        }
        __syncthreads();

        seq_idx += spg;
    }

    states[ bid * tpb + tid ] = local_state;
}

/**
 *
 * Hash-less
 */
template < class StateType, class AlleleSpaceType, class RealType, class IntType >
__global__ void crossover_kernel( StateType * states
                                , AlleleSpaceType * alleles
                                , device_free_space< IntType, unordered_tag > * free_space
                                , poisson_cdf< RealType, 32 > * pois
                                , device_sequence_space< IntType > * sequences
                                , clotho::utility::algo_version< 3 > * v ) {

    typedef StateType                                   state_type;
    typedef AlleleSpaceType                             allele_space_type;
    typedef typename allele_space_type::real_type       real_type;

    typedef device_sequence_space< IntType >            sequence_space_type;
    typedef typename sequence_space_type::int_type      int_type;
    typedef xover_config< unordered_tag, 3 >            xover_type;

    typedef poisson_cdf< RealType, 32 > poisson_type;

    assert( blockDim.y <= xover_type::MAX_WARPS && blockDim.x == 32);   // 8 or fewer full warps

    unsigned int  nAlleles = alleles->capacity;
    if( nAlleles == 0 ) { return; }

    assert( (nAlleles & 31) == 0 ); // nAlleles == 32 * m

    const unsigned int MAX_EVENTS_PER_WARP = 128;    // maximum number of recombination events per sequence
    __shared__ real_type rand_pool[ MAX_EVENTS_PER_WARP * xover_type::MAX_WARPS ];  // 16 warps/block (arbitrary number); 512 random numbers
    __shared__ real_type        s_pois_cdf[ poisson_type::MAX_K ];                      // 4 * 32 == 128

    unsigned int max_k = pois->max_k;

    if( threadIdx.y == 0 ) {    // use first warp to read poisson CDF
        s_pois_cdf[ threadIdx.x ] = pois->_cdf[threadIdx.x];
    }
    __syncthreads();

    unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    unsigned int tpb = (blockDim.x * blockDim.y);
    unsigned int bpg = (gridDim.x * gridDim.y);

    unsigned int spg = bpg * blockDim.y;
    
    int_type seq_width = sequences->seq_width;
    int_type nSequences = sequences->seq_count;    

    int_type    * seqs = sequences->sequences;
    real_type   * allele_list = alleles->locations;

    unsigned int seq_idx = bid * blockDim.y + threadIdx.y;
    state_type local_state = states[ seq_idx * 32 + threadIdx.x ]; // every thread/warp uses the SAME random state

    unsigned int max_seqs = nSequences / blockDim.y;
    max_seqs += (( nSequences % blockDim.y) ? 1 : 0);
    max_seqs *= blockDim.y;

    while( seq_idx < max_seqs ) {
        real_type x = curand_uniform( &local_state );
        unsigned int rand = _find_poisson_maxk32( s_pois_cdf, x, max_k );
        __syncthreads();

        for(unsigned int i = 1; i < 32; i <<= 1 ) {
            unsigned int t = __shfl_up( rand, i );
            rand += ((threadIdx.x >= i) * t);
        }

        unsigned int max_events = __shfl( rand, 31 );
        // fill random pool
        //
        if( max_events >= MAX_EVENTS_PER_WARP ) {
            if( threadIdx.x == 0 ) {
                printf( "Too many events to generate: %d\n", max_events );
            }
            assert( max_events < MAX_EVENTS_PER_WARP );
        }
        __syncthreads();

        rand_pool[ tid ] = curand_uniform( &local_state );
        rand_pool[ tid + tpb ] = curand_uniform( &local_state );
        rand_pool[ tid + 2 * tpb ] = curand_uniform( &local_state );
        rand_pool[ tid + 3 * tpb ] = curand_uniform( &local_state );
        __syncthreads();

        unsigned int seq_offset = seq_idx * seq_width;
        unsigned int a_idx = threadIdx.x;
        while( a_idx < nAlleles ) {
            real_type loc = allele_list[ a_idx ];

            unsigned int s = 0, mask = 0;
            while( s < max_events ) {   // warps within a block diverge. threads within a warp do not
                real_type y = rand_pool[ s++ ]; // every thread within a warp read/loads the same value (sequence offset)
                mask += ( y < loc );
            }
            __syncthreads();

            mask = ((mask & 1) << threadIdx.x);

            for( unsigned int i = 1; i < 32; i <<= 1 ) {
                unsigned int e = __shfl_up( mask, i );
                mask |= ((threadIdx.x >= i) * e);
            }

            if( threadIdx.x == 31 && seq_idx < nSequences ) {
                seqs[ seq_offset ] = mask;
            }
            __syncthreads();
            a_idx += 32;
            ++seq_offset;
        }
        seq_idx += spg;
    }

    seq_idx = bid * blockDim.y + threadIdx.y;  // reset seq_idx
    states[ seq_idx * 32 + threadIdx.x ] = local_state;
}

template < class StateType, class AlleleSpaceType, class RealType, class IntType >
__global__ void crossover_kernel( StateType * states
                                , AlleleSpaceType * alleles
                                , device_free_space< IntType, unordered_tag > * free_space
                                , poisson_cdf< RealType, 32 > * pois
                                , device_sequence_space< IntType > * sequences
                                , clotho::utility::algo_version< 4 > * v ) {

    typedef StateType                                   state_type;
    typedef AlleleSpaceType                             allele_space_type;
    typedef typename allele_space_type::real_type       real_type;

    typedef device_sequence_space< IntType >            sequence_space_type;
    typedef typename sequence_space_type::int_type      int_type;

    typedef poisson_cdf< RealType, 32 >                 poisson_type;
    typedef xover_config< unordered_tag, 4 >            xover_type;

    assert( blockDim.x == 32 && blockDim.y <= xover_type::MAX_WARPS );

    int_type bid = blockIdx.y * gridDim.x + blockIdx.x;
    int_type tid = threadIdx.y * blockDim.x + threadIdx.x;

    unsigned int tpb = (blockDim.x * blockDim.y);
    assert( (tpb % allele_space_type::ALIGNMENT_SIZE) == 0 );

    unsigned int bpg = (gridDim.x * gridDim.y);

    unsigned int i;
    int_type sequence_width = sequences->seq_width;
    int_type nSequences = sequences->seq_count;    

    int_type    * seqs = sequences->sequences;

    //int_type  cap = sequences->capacity;

    real_type       * allele_list = alleles->locations;
    unsigned int    nAlleles = alleles->capacity;

    assert( (nAlleles % allele_space_type::ALIGNMENT_SIZE) == 0 );

    if( nAlleles == 0 ) { return; }

    __shared__ real_type        s_pois_cdf[ poisson_type::MAX_K ];
    __shared__ real_type        rand_pool[ allele_space_type::ALIGNMENT_SIZE ];
    __shared__ unsigned int     event_hash[ allele_space_type::ALIGNMENT_SIZE];

    unsigned int max_k = pois->max_k;
    if( tid < poisson_type::MAX_K ) {
        s_pois_cdf[ tid ] = pois->_cdf[tid];
    }
    __syncthreads();

    state_type local_state = states[ bid * tpb + tid ];

    unsigned int nonzero_warp = (threadIdx.y != 0);
    unsigned int nonzero_thread = (tid != 0);
    int_type seq_idx = bid;
    while( seq_idx < nSequences ) {
        real_type x = curand_uniform( &local_state );   // x in (0, 1]
        rand_pool[ tid ] = ((x >= 1.0) ? 0.0 : x);  // wrap around x to be in [0, 1); this way all event hash bins follow [,) pattern

        x = curand_uniform( &local_state );

        int_type rand = _find_poisson_maxk32( s_pois_cdf, x, max_k );
        __syncthreads();

        for( i = 1; i < 32; i <<= 1 ) {
            unsigned int r = __shfl_up( rand, i );
            rand += ( (threadIdx.x >= i ) * r );
        }

        if( threadIdx.x == 31 ) {
            event_hash[ threadIdx.y ] = rand;
        }
        __syncthreads();

        unsigned int _sum = event_hash[ threadIdx.x ];
        _sum *= (threadIdx.x < blockDim.y);
        __syncthreads();

        for( i = 1; i < 32; i <<= 1 ) {
            unsigned int s = __shfl_up( _sum, i );
            _sum += (( threadIdx.x >= i ) * s);
        }

        unsigned int s = __shfl( _sum, 31 );
//        assert( max_events < allele_space_type::ALIGNMENT_SIZE );
//
        if( s == 0 ) { // true for all threads in block assuming 1 block per sequence
            // if there are no events for this sequence, then simply clear the memory
            unsigned int seq_start = seq_idx * sequence_width + tid;
            unsigned int seq_end = (seq_idx + 1) * sequence_width;
            while( seq_start < seq_end ) {
                seqs[ seq_start ] = 0;
                seq_start += tpb;
            }
            __syncthreads();
        } else {

            s = __shfl( _sum, threadIdx.y - nonzero_warp);
            s *= nonzero_warp;
            __syncthreads();

            rand += s;
            event_hash[tid] = rand;
            __syncthreads();

            i = event_hash[ tid - nonzero_thread];    // minimum event index
            i *= nonzero_thread;
            __syncthreads();

            // BEGIN divergent code
            while (i < rand) {
                x = rand_pool[ i ];
                rand_pool[ i++ ] = (((real_type) tid) + x) / ((real_type) allele_space_type::ALIGNMENT_SIZE);
            }
            __syncthreads();
            // END divergent code

            unsigned int seq_offset = seq_idx * sequence_width + threadIdx.y;   // every thread in a warp has the same block offset
            i = tid;
            while( i < nAlleles ) {
                x = allele_list[ i ];

                rand = (unsigned int) ( x * ((real_type)allele_space_type::ALIGNMENT_SIZE));

                unsigned int nonzero_bin = (rand != 0);  // _keep == 0 -> rand == 0 -> e_min == 0; _keep == 1 -> rand != 0 -> e_min == event_hash[ rand - 1]

                // each thread reads hash (and random pool) relative to their 
                // local allele (x)
                // this code will result in bank conflicts
                //
                // initial performance results suggest that this is
                // an acceptable overhead as overall runtime of simulation
                // loop is minimized when this algorithm is used
                unsigned int e_max = event_hash[ rand ];
                unsigned int e_min = event_hash[ rand - nonzero_bin ];
                e_min *= nonzero_bin;

                int_type cmask = e_min;

                // BEGIN divergent code
                while( e_min < e_max ) {
                    real_type y = rand_pool[ e_min++ ];
                    cmask += (y < x);
                }
                __syncthreads();
                // END divergent code

                cmask = ((cmask & 1) << threadIdx.x);
                
                for( unsigned int j = 1; j < 32; j <<= 1 ) {
                    int_type _c = __shfl_up( cmask, j);
                    cmask |= ((threadIdx.x >= j) * _c);
                }
                if( threadIdx.x == 31 ) {
//                    assert( seq_offset < cap );
                    seqs[ seq_offset ] = cmask;
                }
                __syncthreads();

                i += tpb;
                seq_offset += blockDim.y;   // wpb
            }
            __syncthreads();

        }
        seq_idx += bpg;
    }

    states[ bid * tpb + tid ] = local_state;
}

#endif  // CROSSOVER_KERNEL_UNORDERED_IMPL_HPP_
