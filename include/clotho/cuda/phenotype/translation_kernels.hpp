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
#ifndef TRANSLATE_KERNELS_HPP_
#define TRANSLATE_KERNELS_HPP_

#include <cuda.h>

#include "clotho/cuda/data_spaces/allele_space/device_allele_space.hpp"
#include "clotho/cuda/data_spaces/sequence_space/device_sequence_space.hpp"
//#include "clotho/cuda/data_spaces/phenotype_space/device_phenotype_space.hpp"
#include "clotho/cuda/data_spaces/basic_data_space.hpp"

#include "clotho/utility/algorithm_version.hpp"

/**
 * 
 * Combine and Translate genetic sequence to phenotype
 *
 * 1 block == 1 sequence
 */
template < class RealType, class IntType >
__global__ void _translate( device_weighted_allele_space< RealType > * alleles
                            , device_sequence_space< IntType > * sequences
                            , basic_data_space< RealType > * phenos 
                            , clotho::utility::algo_version< 1 > * v
                            , unsigned int seq_start
                            , unsigned int seq_end) {

    assert( blockDim.x == 32 );     // assert block is warp aligned

    typedef typename device_allele_space< RealType >::real_type     effect_type;
    typedef typename basic_data_space< RealType >::value_type  real_type;
    typedef typename device_sequence_space< IntType >::int_type     int_type;
    typedef device_weighted_allele_space< RealType >    allele_space_type;
    typedef device_sequence_space< IntType >            sequence_space_type;

    unsigned int allele_count = alleles->capacity;
    unsigned int s_width = sequences->seq_width;

    // validation check
    assert( allele_count == s_width * sequence_space_type::OBJECTS_PER_INT );

    int_type    * seqs = sequences->sequences;
    effect_type * effects = alleles->weights;

    __shared__ real_type s_pheno[ 32 ];

    const int_type thread_mask = ((int_type) 1 << threadIdx.x);

    seq_start += blockIdx.y*gridDim.x + blockIdx.x;
    assert( seq_end <= sequences->seq_count );
    while( seq_start < seq_end ) {

        if( threadIdx.y == 0 ) {    // use first warp
            // clear buffer
            // necessary for when blockDim.y != 32
            s_pheno[threadIdx.x] = (real_type)0;
        }
        __syncthreads();

        real_type   pheno = 0.;
        unsigned int a_idx = threadIdx.y * blockDim.x + threadIdx.x;
        unsigned int s_idx = seq_start * s_width + threadIdx.y;
        while( a_idx < allele_count ) {
            effect_type eff = effects[a_idx];   // assumes neutral alleles have 0 weight (handled during allele generation)
            int_type    mask = seqs[ s_idx ];

            eff *= (( (mask & thread_mask) != 0) ? (real_type) 1.0 : (real_type) 0.0);
            pheno += eff;

            a_idx += (blockDim.x * blockDim.y);
            s_idx += blockDim.y;
        }
        __syncthreads();

        for( unsigned int i = 1; i < 32; i <<= 1 ) {
            real_type p = __shfl_up( pheno, i );
            pheno += ((real_type)(threadIdx.x >= i)) * p;
        }

        if( blockDim.y > 1 ) {  // true/false for all threads in block
            if( threadIdx.x == 31 ) {   // write each warp's phenotype partial sum to shared memory
                s_pheno[threadIdx.y] = pheno;
            }
            __syncthreads();

            // Have all warps read the partial phenotype sums
            //  and compute final phenotype sum
            // All warps used to 'avoid' thread divergence.
            real_type tmp_pheno = s_pheno[threadIdx.x];
            for(unsigned int i = 1; i < 32; i <<= 1 ) {
                real_type p = __shfl_up( tmp_pheno, i );
                pheno += ((real_type)(threadIdx.x >= i )) * p;
            }
        }
        __syncthreads();

        if( threadIdx.y == 0 && threadIdx.x == 31 ) {
            phenos->data[ seq_start ] = pheno;
        }
        __syncthreads();

        seq_start += gridDim.x * gridDim.y;
    }
}

/**
 *
 * Translate genetic sequence to phenotype
 * 1 warp == 1 sequence
 */
template < class RealType, class IntType >
__global__ void _translate( device_weighted_allele_space< RealType > * alleles
                            , device_sequence_space< IntType > * sequences
                            , basic_data_space< RealType > * phenos 
                            , clotho::utility::algo_version< 2 > * v
                            , unsigned int seq_start
                            , unsigned int seq_end) {

    assert( blockDim.x == 32 ); // assert that threads are indexed by warps

    assert( (seq_end - seq_start) % blockDim.y == 0 );  // assert no trailing sequences in looping

    typedef device_weighted_allele_space< RealType >    allele_space_type;
    typedef device_sequence_space< IntType >            sequence_space_type;

    typedef typename device_allele_space< RealType >::real_type     effect_type;
    typedef typename basic_data_space< RealType >::value_type       real_type;
    typedef typename device_sequence_space< IntType >::int_type     int_type;

    unsigned int s_width = sequences->seq_width;

    assert( s_width * sequence_space_type::OBJECTS_PER_INT == alleles->capacity );

    int_type    * seqs = sequences->sequences;
    effect_type * effect_sizes = alleles->weights;

    assert( phenos->capacity == sequences->seq_count );

    const int_type thread_mask = ((int_type) 1 << threadIdx.x);

    seq_start +=  (blockIdx.y * gridDim.x + blockIdx.x) * blockDim.y + threadIdx.y;
    while( seq_start < seq_end ) {

        real_type   pheno = 0.;
        unsigned int start = seq_start * s_width;

        for( unsigned int i = 0, j = threadIdx.x; i < s_width; ++i, j += blockDim.x ) {   // all warps process the same number of blocks
            effect_type eff = effect_sizes[ j ];

            int_type    mask = seqs[ start + i ];
            eff *= (( (mask & thread_mask) != 0) ? (real_type) 1.0 : (real_type) 0.0);
            pheno += eff;
        }
        __syncthreads();

        // scan right the phenotype per warp (sequence)
        for( unsigned int i = 1; i < 32; i <<= 1 ) {
            real_type p = __shfl_up( pheno, i );
            pheno += ((real_type)(threadIdx.x >= i)) * p;
        }

        if( threadIdx.x == 31  ) {
            phenos->data[ seq_start ] = pheno;
        }
        __syncthreads();

        seq_start += blockDim.y * gridDim.x * gridDim.y;
    }
}

/**
 *
 * Translate genetic sequence to phenotype
 * 1 thread == 1 sequence
 */
template < class RealType, class IntType >
__global__ void _translate( device_weighted_allele_space< RealType > * alleles
                            , device_sequence_space< IntType > * sequences
                            , basic_data_space< RealType > * phenos 
                            , clotho::utility::algo_version< 3 > * v
                            , unsigned int seq_start
                            , unsigned int seq_end ) {

    assert( blockDim.x == 32 );

    typedef device_weighted_allele_space< RealType >    allele_space_type;
    typedef device_sequence_space< IntType >            sequence_space_type;

    typedef typename device_allele_space< RealType >::real_type     effect_type;
    typedef typename basic_data_space< RealType >::value_type       real_type;
    typedef typename device_sequence_space< IntType >::int_type     int_type;

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int sequences_per_block = blockDim.x * blockDim.y;

    unsigned int s_width = sequences->seq_width;

    int_type    * seqs = sequences->sequences;
    effect_type * effect_sizes = alleles->weights;

    unsigned int padded_end = seq_end + (sequences_per_block - ((seq_end - seq_start) % sequences_per_block));

    assert( (padded_end - seq_start) % sequences_per_block == 0 );

    __shared__ real_type eff_size[ 32 ];

    seq_start +=  (blockIdx.y * gridDim.x + blockIdx.x) * sequences_per_block + tid ;
    __syncthreads();

    while( seq_start < padded_end ) {
        bool write_result = (seq_start < seq_end);

        // the last sequencs should be broadcast to all padded threads.
        if( seq_start >= seq_end ) {
            seq_start = seq_end - 1;
        }
        __syncthreads();

        real_type pheno = 0.0;
        for( unsigned int i = 0; i < s_width; ++i ) {
            // use the first warp to fill the allele intermediate layer
            if( threadIdx.y == 0 ) {
                eff_size[ threadIdx.x ] = effect_sizes[ i * blockDim.x + threadIdx.x ];
            }
            __syncthreads();

            int_type mask = seqs[ seq_start * s_width + i ];
            for( unsigned int j = 0; j < 32; ++j, mask >>= 1 ) {
                real_type eff = eff_size[ j ];
                eff *= (((mask & 1) != 0) ? (real_type) 1.0 : (real_type) 0.0);
                pheno += eff;
            }

            __syncthreads();
        }

        if( write_result ) {
            phenos->data[ seq_start ] = pheno;
        }
        __syncthreads();

        seq_start += (gridDim.x * gridDim.y) * sequences_per_block;
    }

}

#endif  // TRANSLATE_KERNELS_HPP_
