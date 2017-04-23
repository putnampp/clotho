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
                            , clotho::utility::algo_version< 1 > * v) {

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

    unsigned int seq_idx = blockIdx.y*gridDim.x + blockIdx.x;
    const unsigned int seq_count = sequences->seq_count;
    while( seq_idx < seq_count ) {

        if( threadIdx.y == 0 ) {    // use first warp
            // clear buffer
            // necessary for when blockDim.y != 32
            s_pheno[threadIdx.x] = (real_type)0;
        }
        __syncthreads();

        real_type   pheno = 0.;
        unsigned int a_idx = threadIdx.y * blockDim.x + threadIdx.x;
        unsigned int s_idx = seq_idx * s_width + threadIdx.y;
        while( a_idx < allele_count ) {
            effect_type eff = effects[a_idx];   // assumes neutral alleles have 0 weight (handled during allele generation)
            int_type    mask = seqs[ s_idx ];

            pheno += ((real_type)((mask >> threadIdx.x) & 1)) * eff;

            a_idx += (blockDim.x * blockDim.y);
            s_idx += blockDim.y;
        }
        __syncthreads();

        for( unsigned int i = 1; i < 32; i <<= 1 ) {
            real_type p = __shfl_up( pheno, i );
            pheno += ((real_type)(threadIdx.x >= i)) * p;
        }

        if( warp_per_block > 1 ) {  // true/false for all threads in block
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
            phenos->data[ seq_idx ] = pheno;
        }
        __syncthreads();

        seq_idx += gridDim.x * gridDim.y;
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
                            , clotho::utility::algo_version< 2 > * v) {

    typedef device_weighted_allele_space< RealType >    allele_space_type;
    typedef device_sequence_space< IntType >            sequence_space_type;

    typedef typename device_allele_space< RealType >::real_type     effect_type;
    typedef typename basic_data_space< RealType >::value_type  real_type;
    typedef typename device_sequence_space< IntType >::int_type     int_type;

    unsigned int allele_count = alleles->capacity;
    unsigned int s_width = sequences->seq_width;

    assert( allele_count == s_width * sequence_space_type::OBJECTS_PER_INT );

    assert( ((blockDim.x * blockDim.y) & 31) == 0 );  // all warps are full

    unsigned int wpb = ((blockDim.x * blockDim.y) >> 5);  // warp (==phenotype==individual)/block
    unsigned int bpg = (gridDim.x * gridDim.y); // block/grid
    unsigned int ppg = wpb * bpg;   // phenotypes (==individuals)/grid

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int lane_id = (tid & 31);

    int_type    * seqs = sequences->sequences;
    effect_type * effect_sizes = alleles->weights;

    unsigned int pheno_count = phenos->capacity;
    unsigned int seq_count = sequences->seq_count;

    assert( (seq_count & 1) == 0 ); // even number of sequences
    assert( (2 * pheno_count) == seq_count );   // reduce 2 sequences to 1 phenotype

    unsigned int max_pheno_count = pheno_count / wpb;
    max_pheno_count += ((pheno_count % wpb) ? 1 : 0);
    max_pheno_count *= wpb;

    unsigned int pheno_idx =  (blockIdx.y * gridDim.x + blockIdx.x) * wpb + (tid >> 5);
    while( pheno_idx < max_pheno_count ) {

        real_type   pheno = 0.;
        unsigned int a_idx = lane_id;

        unsigned int end = (2 * pheno_idx + 1) * s_width;
        unsigned int start = end - (( pheno_idx < pheno_count ) ? (s_width) : 0);

        while( start < end ) {
            effect_type eff = effect_sizes[a_idx];

            int_type    mask = seqs[ start ];
            pheno += ((real_type)((mask >> lane_id) & 1)) * eff;

            mask = seqs[ start + s_width ];
            pheno += ((real_type)((mask >> lane_id) & 1)) * eff;

            a_idx += 32;
            ++start;
        }
        __syncthreads();

        for( unsigned int i = 1; i < 32; i <<= 1 ) {
            real_type p = __shfl_up( pheno, i );
            pheno += ((real_type)(lane_id >= i)) * p;
        }

        if( lane_id == 31 && (pheno_idx < pheno_count) ) {
            phenos->data[ pheno_idx ] = pheno;
        }
        __syncthreads();

        pheno_idx += ppg;
    }
}

///**
// *
// * Translate genetic sequence to phenotype
// * 1 warp == 1 sequence
// */
//template < class RealType, class IntType >
//__global__ void _translate( device_weighted_allele_space< RealType > * alleles
//                            , device_sequence_space< IntType > * sequences
//                            , basic_data_space< RealType > * phenos 
//                            , clotho::utility::algo_version< 3 > * v) {
//
//    assert( blockDim.x == 32 ); // assert block dimensionality defined in terms of warps
//
//    while( seq_idx < seq_count ) {
//
//        while( all_idx < all_count ) {
//
//        }
//
//    }
//
//}
#endif  // TRANSLATE_KERNELS_HPP_
