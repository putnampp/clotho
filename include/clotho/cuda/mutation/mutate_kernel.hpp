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
__global__ void mutate( IntType * offspring, unsigned int * allele_index_pool, unsigned int sequence_dist, unsigned int width ) {

    unsigned int start = sequence_dist[ blockDim.x ];
    unsigned int end = sequence_dist[ blockDim.x + 1 ];

    while( start < end ) {
        unsigned int idx = allele_index_pool[ start ];
        IntType mask = ( 1 << (idx % 8 * sizeof(IntType)) );

        unsigned int b_idx = (idx / (8 * sizeof(IntType)));

        b_idx += blockIdx.x * width;

        IntType b = offspring[ b_idx ];

        // assert that not colliding an existing mutation in sequence
        assert( ( b & mask ) == 0 );
        b |= mask;

        offspring[ b_idx ] = b;
        ++start;
    }
}

#endif  // MUTATE_KERNEL_HPP_
