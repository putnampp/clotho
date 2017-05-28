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
#ifndef PHENOTYPE_EVALUATOR_HPP_
#define PHENOTYPE_EVALUATOR_HPP_

/**
 *
 * 1 block per sequence per trait
 *
 */
template < class IntType, class RealType >
__global__ void evaluate_phenotype( IntType * seqs, RealType * trait_weights, RealType * phenos, unsigned int width, unsigned int allele_count ) {

    // gridDim.x == sequence_count
    // gridDim.y == trait_count
    // blockIdx.x == sequence index
    // blockIdx.y == trait index
    // blockDim.x == alleles per sequence block
    // blockDim.y == sequence blocks per batch (warps)
    // threadIdx.x == allele index in sequence block
    // threadIdx.y == offset of sequence block
    unsigned int allele_idx = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int b_idx = threadIdx.y;

    IntType mask = (1 << threadIdx.x);
    RealType local_pheno = 0.0;

    __shared__ RealType sBuffer[ 32 ];

    // clear shared space
    // some what unnecessary when 32 blocks (warps) are in use
    if( threadIdx.y == 0 ) {
        sBuffer[ threadIdx.x ] = 0.0;
    }
    __syncthreads();

    while( allele_idx < allele_count ) {
        RealType w = trait_weights[ blockIdx.y * allele_count + allele_idx ];
        IntType b = seqs[ blockIdx.x * width + b_idx];

        if( (b & m) ) {
            local_pheno += w;
        }

        allele_idx += (blockDim.x * blockDim.y);
        b_idx += blockDim.y;
    }
    __syncthreads();

    for( unsigned int i = 1; i < 32; i <<= 1 ) {
        RealType _p = __shfl_up( local_pheno, i );
        local_pheno += (( threadIdx.x >= i ) ? _p : 0.0);
    }

    if( threadIdx.x == 31 ) {
        sBuffer[ threadIdx.y ] = local_pheno;
    }
    __syncthreads();

    local_pheno = sBuffer[ threadIdx.x ];

    for( unsigned int i = 1; i < 32; i <<= 1 ) {
        RealType _p = __shfl_up( local_pheno, i );
        local_pheno += (( threadIdx.x >= i) ? _p : 0.0);
    }

    if( threadIdx.y == 0 && threadIdx.x == 31 ) {
        phenos[ blockIdx.y * gridDim.x + blockIdx.x ] = local_pheno; 
    }
}

#endif  // PHENOTYPE_EVALUATOR_HPP_
