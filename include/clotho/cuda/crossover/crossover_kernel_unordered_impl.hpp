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

#include "clotho/cuda/data_spaces/allele_space/device_allele_space.hpp"
#include "clotho/cuda/data_spaces/sequence_space/device_sequence_space.hpp"
#include "clotho/cuda/distributions/poisson_distribution.hpp"

#include "clotho/cuda/data_spaces/tags/unordered_tag.hpp"
#include "clotho/cuda/crossover/persist_sequence.hpp"

template < class StateType, class RealType, class IntType >
__global__ void crossover_kernel( StateType * states
                                , device_allele_space< RealType, IntType, unordered_tag > * alleles
                                , poisson_cdf< RealType, 32 > * pois
                                , device_sequence_space< IntType > * sequences ) {

    typedef device_allele_space< RealType, IntType, unordered_tag > allele_space_type;
    typedef RealType   real_type;
    typedef StateType  state_type;
    typedef IntType    int_type;

    typedef poisson_cdf< RealType, 32 > poisson_type;

    int_type tid = threadIdx.y * blockDim.x + threadIdx.x;
    int_type lane_id = (tid & 31);
    int_type warp_id = (tid >> 5);

    int_type i;
    int_type sequence_width = sequences->seq_width;
    int_type nSequences = sequences->seq_count;    

    int_type    * seqs = sequences->sequences;

    real_type   * allele_list = alleles->locations;
    int_type    nAlleles = alleles->size;

    if( nAlleles == 0 ) { return; }

    __shared__ real_type        s_pois_cdf[ poisson_type::MAX_K ];
    __shared__ real_type        rand_pool[ allele_space_type::ALIGNMENT_SIZE ];
    __shared__ unsigned int     event_hash[ allele_space_type::ALIGNMENT_SIZE + 1];

    event_hash[ tid ] = 0;

    unsigned int max_k = pois->max_k;
    if( tid < poisson_type::MAX_K ) {
        s_pois_cdf[ tid ] = pois->_cdf[tid];
    }
    __syncthreads();

    state_type local_state = states[ blockIdx.x * blockDim.x *blockDim.y + tid ];

    int_type seq_idx = blockIdx.x;
    while( seq_idx < nSequences ) {

        int_type * seq = seqs + (seq_idx * sequence_width);

        real_type x = curand_uniform( &local_state );
        rand_pool[ tid ] = curand_uniform( &local_state );

        int_type rand = _find_poisson_maxk32( s_pois_cdf, x, max_k );
        __syncthreads();

        for( i = 1; i < 32; i <<= 1 ) {
            unsigned int r = __shfl_up( rand, i );
            rand += ( (lane_id >= i ) * r );
        }

        if( lane_id == 31 ) {
            event_hash[ 32 + warp_id ] = rand;
        }
        __syncthreads();

        unsigned int _sum = event_hash[ 32 + lane_id ];
        for( i = 1; i < (allele_space_type::ALIGNMENT_SIZE >> 5); i <<= 1 ) {
            unsigned int s = __shfl_up( _sum, i );
            _sum += (( lane_id >= i ) * s);
        }

        unsigned int s = __shfl( _sum, warp_id - 1);
        rand += (( warp_id != 0 ) * s);
        __syncthreads();

        event_hash[ tid + 1 ] = rand;
        __syncthreads();

        i = event_hash[tid];    // minimum event index
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


        i = tid;
        while( i < nAlleles ) {
            x = allele_list[ i ];

            rand = (unsigned int) ( x * ((real_type)allele_space_type::ALIGNMENT_SIZE));

            unsigned int e_min = event_hash[ rand++ ];
            _sum = event_hash[ rand ];

            int_type cmask = e_min;

            // BEGIN divergent code
            while( e_min < _sum ) {
                accum = rand_pool[ e_min++ ];
                cmask += ( x > accum);
            }
            __syncthreads();
            // END divergent code

            cmask = ((cmask & 1) * (1 << lane_id));

            persist_mask_unrolled( cmask, warp_id, lane_id, seq );
            __syncthreads();

            i += allele_space_type::ALIGNMENT_SIZE;
            seq += (allele_space_type::ALIGNMENT_SIZE >> 5);
        }
        __syncthreads();


        seq_idx += (gridDim.x * gridDim.y);
    }

    states[ blockIdx.x * blockDim.x *blockDim.y + tid ] = local_state;
}

#endif  // CROSSOVER_KERNEL_UNORDERED_IMPL_HPP_
