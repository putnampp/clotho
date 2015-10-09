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
#ifndef CROSSOVER_KERNEL_UNIT_ORDERED_IMPL_HPP_
#define CROSSOVER_KERNEL_UNIT_ORDERED_IMPL_HPP_

#include "clotho/cuda/data_spaces/allele_space/device_allele_space.hpp"
#include "clotho/cuda/data_spaces/sequence_space/device_sequence_space.hpp"
#include "clotho/cuda/data_spaces/free_space/device_free_space.hpp"
#include "clotho/cuda/distributions/poisson_distribution.hpp"

#include "clotho/cuda/data_spaces/tags/unit_ordered_tag.hpp"
#include "clotho/utility/algorithm_version.hpp"

template < class StateType, class AlleleSpaceType, class RealType, class IntType, unsigned char V >
__global__ void crossover_kernel( StateType * states
                                , AlleleSpaceType * alleles
                                , device_free_space< IntType, unit_ordered_tag< IntType > > * free_space
                                , poisson_cdf< RealType, 32 > * pois
                                , device_sequence_space< IntType > * offspring_seqs
                                , clotho::utility::algo_version< V > * v ) {

    typedef StateType                                       state_type;

    typedef poisson_cdf< RealType, 32 >                     poisson_type;
    typedef device_sequence_space< IntType >        sequence_space_type;
    typedef typename sequence_space_type::int_type          sequence_int_type;

    typedef AlleleSpaceType                                 allele_space_type;
    typedef typename allele_space_type::real_type           real_type;

    unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;
    unsigned int offspring_count = offspring_seqs->seq_count;

    if( bid >= offspring_count ) return;

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int warp_id = (tid >> 5 );
    unsigned int lane_id = (tid & 31 );

    unsigned int tpb = (blockDim.x * blockDim.y);
    unsigned int wpb = (tpb >> 5);
    unsigned int bpg = (gridDim.x * gridDim.y);
    unsigned int off_width = offspring_seqs->seq_width;
    unsigned int allele_cap = alleles->capacity;    // capacity is a multiple block dimensions (thread count)

    state_type local_state = states[ bid * tpb + tid ];

    __shared__ unsigned int s_event_hash[ 32 ];
    __shared__ real_type    s_rand_pool[ 1024 ];
    __shared__ real_type    s_pois_cdf[ poisson_type::MAX_K ];

    unsigned int max_k = pois->max_k;
    if( tid < poisson_type::MAX_K ) {
        s_pois_cdf[ tid ] = pois->_cdf[ tid ];
    }
    __syncthreads();

    real_type * all_loc = alleles->locations;
    sequence_int_type * off_seqs = offspring_seqs->sequences;

    while( bid < offspring_count ) {    // true for all threads
        // every thread generates a random number
        real_type r = curand_uniform( &local_state );

        unsigned int e = _find_poisson_maxk32( s_pois_cdf, r, max_k );
        __syncthreads();

        if( warp_id == 0 ) {    // only save first warps random numbers
            s_event_hash[ lane_id ] = e;
        }
        __syncthreads();

        e = s_event_hash[ lane_id ];    // broadcast event count by lane to all warps
        __syncthreads();
        
        unsigned int psum = e;          // preserve the lane's event count for partial summing
        
        // determine the maximum number of events of all lanes within a warp
        for( unsigned int i = 1; i < 32; i <<= 1 ) {
            unsigned int _e = __shfl_up( e, i );
            e = ((_e > e) ? _e : e );
        }
        unsigned int max = __shfl( e, 31 ); // all warps will share the same max value

        // compute the prefix sum of the event counts with each warp
        for( unsigned int i = 1; i < 32; i <<= 1 ) {
            unsigned int p = __shfl_up( psum, i );
            psum += (( lane_id >= i ) * p );
        }

        unsigned int base_count = __shfl_up( psum, 1 );
        if( lane_id == 0 ) base_count = 0;
        __syncthreads();

        for( unsigned int i = warp_id; i < max; i += wpb ) {    // allows warp divergence
            // over fills the random event pool, but all threads within a warp are active
            r = curand_uniform( &local_state );

            s_rand_pool[ i * 32 + lane_id ] = (((real_type) lane_id + r) / (real_type) 32);
        }
        __syncthreads();

        unsigned int o_idx = bid * off_width + warp_id;
        unsigned int a_id = tid;

        while( a_id < allele_cap ) {    // true for all threads
            real_type a = all_loc[ a_id ];  // every thread reads a location just in case

            e = (psum - base_count);
            unsigned int x = base_count;
            unsigned int use_a = ( a_id < allele_cap );

            for( unsigned int i = 0; i < max; ++i ) {
                // similar to the over filling of the random pool
                // threads where (e < max) will over compute
                // that is, when e <= i < max, x += (0 * ?);
                // in this scenario reading/loading y and performing (y < a)
                // however, allowing this avoids thread divergence

                real_type y = s_rand_pool[ i * 32 + lane_id ];
                x += ((i < e) & use_a & (y < a));   // x += 1 or 0
            }
            __syncthreads();

            x = ((x & 1) << lane_id);

            // shuffle up mask
            for( unsigned int i = 1; i < 32; i <<= 1 ) {
                unsigned int _x = __shfl_up( x, i );
                x |= ((lane_id >= i ) * _x );
            }

            if( lane_id == 31 ) {
                off_seqs[ o_idx ] = x;
            }
            __syncthreads();

            o_idx += wpb;
            a_id += tpb;
        }
        bid += bpg;
    }

    bid = blockIdx.y * gridDim.x + blockIdx.x;
    states[ bid * tpb + tid ] = local_state;
}

template < class StateType, class AlleleSpaceType, class RealType, class IntType >
__global__ void crossover_kernel( StateType * states
                                , AlleleSpaceType * alleles
                                , device_free_space< IntType, unit_ordered_tag< IntType > > * free_space
                                , poisson_cdf< RealType, 32 > * pois
                                , device_sequence_space< IntType > * offspring_seqs
                                , clotho::utility::algo_version< 2 > * v ) {

    typedef StateType                                       state_type;

    typedef poisson_cdf< RealType, 32 >                     poisson_type;
    typedef device_sequence_space< IntType >                sequence_space_type;
    typedef typename sequence_space_type::int_type          sequence_int_type;

    typedef AlleleSpaceType                                 allele_space_type;
    typedef typename allele_space_type::real_type           real_type;

    assert( blockDim.x == 32 ); // every warp is full

    real_type           * allele_loc = alleles->locations;
    unsigned int        allele_count = alleles->capacity;
    assert( (allele_count & 31) == 0);  // allele_count == 32*n

    unsigned int width = offspring_seqs->seq_width;
    if( width == 0 ) return;
    assert( width * sequence_space_type::OBJECTS_PER_INT == allele_count);

//    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    const unsigned int MAX_EVENTS_PER_LANE = poisson_type::MAX_K;   // 32
    real_type l_rand_pool[ MAX_EVENTS_PER_LANE ];
    __shared__ real_type s_pois_cdf[ poisson_type::MAX_K ];

    unsigned int max_k = pois->max_k;
    if( threadIdx.y == 0 ) {
        s_pois_cdf[ threadIdx.x ] = pois->_cdf[ threadIdx.x ];
    }
    __syncthreads();

    const unsigned int spg = (blockDim.y * (gridDim.x * gridDim.y)); // sequences per grid

    sequence_int_type   * seqs = offspring_seqs->sequences;
    unsigned int seq_count = offspring_seqs->seq_count;

    unsigned int state_idx = (blockIdx.y * gridDim.x + blockIdx.x) * (blockDim.x*blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x);
    state_type local_state = states[ state_idx ];   // poor global memory access; seems to perform the equivalent memcpy operation

    unsigned int max_seq_id = seq_count / blockDim.y;
    max_seq_id += ((seq_count % blockDim.y) ? 1 : 0);
    max_seq_id *= blockDim.y;

    unsigned int seq_id = (blockIdx.y * gridDim.x + blockIdx.x) * blockDim.y + threadIdx.y;
    while( seq_id < max_seq_id ) {  // allows blocks to terminate early
        real_type r = curand_uniform( &local_state );

        unsigned int hi = _find_poisson_maxk32( s_pois_cdf, r, max_k );

        unsigned int m = hi;
        for( unsigned int i = 1; i < 32; i <<= 1 ) {
            unsigned int _m = __shfl_up( m, i);
            m = ((_m > m)? _m : m);
        }

        unsigned int max_events = __shfl( m, 31);
        // assert( max_events < MAX_EVENTS_PER_LANE );  // should terminate every thread in a warp; uncertain whether this will terminate the kernel (if a single thread fails an assertion, does the entire kernel execution fail?)

        for( unsigned int i = 1; i < 32; i <<= 1 ) {
            unsigned int e = __shfl_up( hi, i );
            hi += ((threadIdx.x >= i) * e);
        }

        unsigned int lo = __shfl_up( hi, 1 );
        lo *= (threadIdx.x != 0);   // make sure lo = 0 when lane_id == 0

        unsigned int s = 0, ecount = lo;
        while( s < max_events ) {   // all threads in a warp execute same number of steps
            // warps within a block will diverge
            real_type x = curand_uniform( &local_state );
            x = ((ecount++ < hi) ? x : 1.0);
            l_rand_pool[ s++ ] = ((real_type)threadIdx.x + x)/((real_type)32);
        }
        __syncthreads();

        unsigned int idx = seq_id * width;
        unsigned int a_id = threadIdx.x;
        while( a_id < allele_count ) {  // all warps execute the same loop
            real_type x = allele_loc[ a_id ];

            ecount = lo;
            s = 0;
            while( s < max_events ) {   // threads within a warp will not diverge
                real_type y = l_rand_pool[ s++ ];   // rand_pool is in each threads local memory
                ecount += ( y < x );
            }
            __syncthreads();

            ecount = ((ecount & 1) << threadIdx.x);

            for(unsigned int i = 1; i < 32; i <<= 1 ) {
                unsigned int tmp = __shfl_up( ecount, i );
                ecount |= ((threadIdx.x >= i) * tmp);
            }

            if( threadIdx.x == 31 && seq_id < seq_count ) {
                seqs[ idx ] = ecount;
            }
            __syncthreads();
            a_id += 32;
            ++idx;
        }
        seq_id += spg;
    }
    states[ state_idx ] = local_state;
}

#endif  // CROSSOVER_KERNEL_UNIT_ORDERED_IMPL_HPP_
