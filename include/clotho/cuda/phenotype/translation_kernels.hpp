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
#include "clotho/cuda/data_spaces/phenotype_space/device_phenotype_space.hpp"

template < class RealType, class IntType >
__global__ void _translate( device_weighted_allele_space< RealType > * alleles
                            , device_sequence_space< IntType > * sequences
                            , device_phenotype_space< RealType > * phenos ) {

    unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;

    unsigned int seq_idx = 2 * bid;
    unsigned int seq_count = sequences->seq_count;

    if( seq_idx >= seq_count ) return;

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    unsigned int warp_id = (tid >> 5);
    unsigned int lane_id = tid & 31;

    unsigned int allele_count = alleles->size;
    unsigned int s_width = sequences->seq_width;

    typedef typename device_allele_space< RealType >::real_type     effect_type;
    typedef typename device_phenotype_space< RealType >::real_type  real_type;
    typedef typename device_sequence_space< IntType >::int_type     int_type;

    int_type    * seqs = sequences->sequences;
    effect_type * effects = alleles->weights;

    unsigned int thread_per_block = (blockDim.x * blockDim.y);
    unsigned int warp_per_block = (thread_per_block >> 5);
    unsigned int seq_per_grid = (gridDim.x * gridDim.y) * 2;

    while( seq_idx < seq_count ) {

        real_type   pheno = 0.;
        unsigned int a_idx = tid;
        unsigned int s_idx = seq_idx * s_width + warp_id;
        while( a_idx < allele_count ) {
            effect_type eff = effects[a_idx];

            int_type    mask = seqs[ s_idx ];
            pheno += ((real_type)((mask >> lane_id) & 1)) * eff;

            mask = seqs[ s_idx + s_width ];
            pheno += ((real_type)((mask >> lane_id) & 1)) * eff;

            a_idx += thread_per_block;
            s_idx += warp_per_block;
        }
        __syncthreads();

        for( unsigned int i = 1; i < 32; i <<= 1 ) {
            real_type p = __shfl_up( pheno, i );
            pheno += ((real_type)(lane_id >= i)) * p;
        }

        if( warp_per_block > 1 ) {  // true/false for all threads in block
            __shared__ real_type s_pheno[ 32 ];
            if( lane_id < 32 ) {
                s_pheno[lane_id] = (real_type)0;
            }
            __syncthreads();
            if( lane_id == 31 ) {
                s_pheno[warp_id] = pheno;
            }
            __syncthreads();

            real_type tmp_pheno = s_pheno[lane_id];
            for(unsigned int i = 1; i < 32; i <<= 1 ) {
                real_type p = __shfl_up( tmp_pheno, i );
                pheno += ((real_type)(lane_id >= i )) * p;
            }
        }
        __syncthreads();

        if( tid == 31 ) {
            phenos->phenotypes[ seq_idx / 2 ] = pheno;
        }
        __syncthreads();

        seq_idx += seq_per_grid;
    }
}

#endif  // TRANSLATE_KERNELS_HPP_
