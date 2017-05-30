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
#ifndef MUTATE_KERNEL_HPP_
#define MUTATE_KERNEL_HPP_

template < class IntType >
__global__ void mutate_kernel( IntType * offspring, unsigned int * allele_index_pool, unsigned int * sequence_dist, unsigned int width ) {

    unsigned int start = sequence_dist[ blockIdx.x ];
    unsigned int end = sequence_dist[ blockIdx.x + 1 ];

    const unsigned int BITS_PER_BLOCK = (8 * sizeof(IntType) );

    while( start < end ) {
        unsigned int idx = allele_index_pool[ start ];
        IntType mask = ( 1 << (idx % BITS_PER_BLOCK ));

        assert( (idx / BITS_PER_BLOCK) < width );
        unsigned int b_idx = blockIdx.x * width + (idx / BITS_PER_BLOCK);

        // read the block
        IntType b = offspring[ b_idx ];

        // assert that not colliding an existing mutation in sequence
        assert( ( b & mask ) == 0 );
        // update the block
        b = (b | mask);

        // write the block 
        offspring[ b_idx ] = b;
        ++start;
    }
}

struct mutate {

    template < class IntType >
    static void execute( IntType * seqs, unsigned int * allele_index_pool, unsigned int * sequence_dist, unsigned int seq_count, unsigned int seq_width ) {
        assert( seqs != NULL );
        assert( allele_index_pool != NULL );
        assert( sequence_dist != NULL );

        dim3 blocks( seq_count, 1, 1 ), threads( 1,1,1 );

//        std::cerr << "Mutate - [" << blocks.x << ", " << blocks.y << "]; [" << threads.x << ", " << threads.y << "]" << std::endl;
        mutate_kernel<<< blocks, threads >>>( seqs, allele_index_pool, sequence_dist, seq_width );
    }
};

#endif  // MUTATE_KERNEL_HPP_
