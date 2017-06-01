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
__global__ void evaluate_phenotype_kernel( IntType * seqs, RealType * trait_weights, RealType * phenos, unsigned int width, unsigned int allele_count, unsigned int trait_width ) {
    assert( blockDim.x == 32 );

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
        // trait_index * allele_count + allele_index
        RealType w = trait_weights[ blockIdx.y * trait_width + allele_idx ];
        // sequence_index * blocks_per_sequence + block_index
        IntType b = seqs[ blockIdx.x * width + b_idx];

        if( (b & mask) != 0) {
            local_pheno += w;
        }

        // step == alleles_per_sequence_block * sequence_blocks
        allele_idx += (blockDim.x * blockDim.y);
        // sequence_blocks
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
        // trait_index * sequence_count + sequence_index
        phenos[ blockIdx.y * gridDim.x + blockIdx.x ] = local_pheno; 
    }
}

template < class IntType, class RealType >
__global__ void evaluate_phenotype_kernel( IntType * seqs, RealType * trait_weights, RealType * phenos, unsigned int seq_count, unsigned int width, unsigned int allele_count, unsigned int trait_width ) {
    assert( blockDim.x == 32 );

    // gridDim.x == sequence_count
    // gridDim.y == trait_count
    // blockIdx.x == sequence index
    // blockIdx.y == trait index
    // blockDim.x == alleles per sequence block
    // blockDim.y == sequence blocks per batch (warps)
    // threadIdx.x == allele index in sequence block
    // threadIdx.y == offset of sequence block
    unsigned int all_idx = threadIdx.x;

    unsigned int seq_idx = blockIdx.x * blockDim.y + threadIdx.y;

    unsigned int b_idx = seq_idx;

    // enables all thread y axis thread groups (warps) to continue to operate
    // despite being outside the bounds of the sequence space
    if( b_idx >= seq_count ) {
        b_idx = seq_count - 1;
    }
    __syncthreads();

    b_idx *= width;

    unsigned int trait_offset = blockIdx.y * allele_count + all_idx;

    IntType mask = (1 << threadIdx.x);
    RealType local_pheno = 0.0;

    unsigned int i = 0;
    while( i  < width ) { // true for all threads in warp and thread block
        // trait_index * allele_count + allele_index
        RealType w = 0.0;

        if( all_idx < allele_count ) {
            w = trait_weights[ trait_offset ];
        }
        __syncthreads();

        // all threads in the same y axis (warp) read the same block
        IntType b = seqs[ b_idx + i ];

        if( (b & mask) != 0) {
            local_pheno += w;
        }

        // step == threads per warp
        all_idx += blockDim.x;
        trait_offset += blockDim.x;

        i += 1;
    }
    __syncthreads();

    for( unsigned int i = 1; i < 32; i <<= 1 ) {
        RealType _p = __shfl_up( local_pheno, i );
        local_pheno += (( threadIdx.x >= i ) ? _p : 0.0);
    }

    if( threadIdx.x == 31 && seq_idx < seq_count ) {
        // trait_index * sequence_count + sequence_index
        phenos[ blockIdx.y * seq_count + seq_idx ] = local_pheno; 
    }
}

#ifndef PHENOTYPE_KERNEL_ALGORITHM
#define PHENOTYPE_KERNEL_ALGORITHM

#define PHENOTYPE_DEFAULT_KERNEL 0
#define PHENOTYPE_WARP_SEQUENCE_KERNEL 1
#endif  // PHENOTYPE_KERNEL_ALGORITHM

struct evaluate_phenotype {

    typedef clotho::utility::algo_version< PHENOTYPE_DEFAULT_KERNEL > default_kernel;
    typedef clotho::utility::algo_version< PHENOTYPE_WARP_SEQUENCE_KERNEL > warp_sequence_kernel;

    template < unsigned char V >
    static void computeDimensions( unsigned int seq_count, unsigned int seq_width, unsigned int trait_count, dim3 & blocks, dim3 & threads, clotho::utility::algo_version< V > * v ) {
        blocks.x = seq_count;
        blocks.y = trait_count;
        if( seq_width < 32 ) {
            threads.y = seq_width;
        }
    }

    static void computeDimensions( unsigned int seq_count, unsigned int seq_width, unsigned int trait_count, dim3 & blocks, dim3 & threads, clotho::utility::algo_version< PHENOTYPE_WARP_SEQUENCE_KERNEL > * v ) {
        blocks.y = trait_count;
        if( seq_count > 32 ) {
            blocks.x = seq_count / 32;
            if( seq_count % 32 ) {
                blocks.x += 1;
            }
            threads.x = 32;
            threads.y = 32;
        } else {
            blocks.x = 1;
            threads.x = 32;
            threads.y = seq_count;
        }
    }

    template < class IntType, class RealType, unsigned char V >
    static void execute( IntType * seqs, RealType * trait_weights, RealType * phenos, unsigned int seq_count, unsigned int seq_width, unsigned int allele_count, unsigned int trait_width, unsigned int trait_count, clotho::utility::algo_version< V > * v ) {
        assert( seqs != NULL );
        assert( trait_weights != NULL );
        assert( phenos != NULL );

        assert( allele_count <= trait_width );
        assert( allele_count <= seq_width * 32 );

        dim3 blocks( 1, 1, 1 ), threads( 32, 32, 1 );
        computeDimensions( seq_count, seq_width, trait_count, blocks, threads, v );

        evaluate_phenotype_kernel<<< blocks, threads >>>( seqs, trait_weights, phenos, seq_width, allele_count, trait_width );
    }

    template < class IntType, class RealType, unsigned char V >
    static void execute( IntType * seqs, RealType * trait_weights, RealType * phenos, unsigned int seq_count, unsigned int seq_width, unsigned int allele_count, unsigned int trait_width, unsigned int trait_count, cudaStream_t & stream, clotho::utility::algo_version< V > * v ) {
        assert( seqs != NULL );
        assert( trait_weights != NULL );
        assert( phenos != NULL );

        assert( allele_count <= trait_width );
        assert( allele_count <= seq_width * 32 );

        dim3 blocks( 1, 1, 1 ), threads( 32, 32, 1 );
        computeDimensions( seq_count, seq_width, trait_count, blocks, threads, v );

        evaluate_phenotype_kernel<<< blocks, threads, 0, stream >>>( seqs, trait_weights, phenos, seq_width, allele_count, trait_width );
    }

    template < class IntType, class RealType >
    static void execute( IntType * seqs, RealType * trait_weights, RealType * phenos, unsigned int seq_count, unsigned int seq_width, unsigned int allele_count, unsigned int trait_width, unsigned int trait_count, clotho::utility::algo_version< PHENOTYPE_WARP_SEQUENCE_KERNEL > * v ) {
        assert( seqs != NULL );
        assert( trait_weights != NULL );
        assert( phenos != NULL );

        assert( allele_count <= trait_width );
        assert( allele_count <= seq_width * 32 );

        dim3 blocks( 1, 1, 1 ), threads( 32, 32, 1 );
        computeDimensions( seq_count, seq_width, trait_count, blocks, threads, v );

        evaluate_phenotype_kernel<<< blocks, threads >>>( seqs, trait_weights, phenos, seq_count, seq_width, allele_count, trait_width );
    }

    template < class IntType, class RealType >
    static void execute( IntType * seqs, RealType * trait_weights, RealType * phenos, unsigned int seq_count, unsigned int seq_width, unsigned int allele_count, unsigned int trait_width, unsigned int trait_count, cudaStream_t & stream, clotho::utility::algo_version< PHENOTYPE_WARP_SEQUENCE_KERNEL > * v ) {
        assert( seqs != NULL );
        assert( trait_weights != NULL );
        assert( phenos != NULL );

        assert( allele_count <= trait_width );
        assert( allele_count <= seq_width * 32 );

        dim3 blocks( 1, 1, 1 ), threads( 32, 32, 1 );
        computeDimensions( seq_count, seq_width, trait_count, blocks, threads, v );

        evaluate_phenotype_kernel<<< blocks, threads, 0, stream >>>( seqs, trait_weights, phenos, seq_count, seq_width, allele_count, trait_width );
    }
};

#endif  // PHENOTYPE_EVALUATOR_HPP_
