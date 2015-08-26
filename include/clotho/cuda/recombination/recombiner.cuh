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
#ifndef RECOMBINER_CUH_
#define RECOMBINER_CUH_


/**
 *  Inline recombiner - assumes crossover pattern stored in child_seqs and
 *  replaces crossover pattern with recombined parental sequence
 *  
 *  Execution Configuration:
 *      - 1 block/child sequence
 *      - 1 parent_index/child sequence (assumes that parental sequences are adjacently indexed)
 *      - blockIdx.x == child_sequence offset
 *      
 */
template < class IntType >
__global__ void diploid_recombiner( IntType * parent_seqs
                            , IntType * parent_index
                            , IntType * child_seqs
                            , IntType parent_sequence_width
                            , IntType child_sequence_width ) {
    IntType tid = threadIdx.y * blockDim.x + threadIdx.x;

    IntType p_idx = parent_idx[ blockIdx.x ];
    p_idx *= 2;
    __syncthreads();

    IntType p = 0, q = 0, c = 0;
    if( tid < parent_sequence_width ) {
        p = parent_seqs[ p_idx * parent_sequence_width + tid ];
        q = parent_seqs[ (p_idx + 1) * parent_sequence_width + tid ];
        c = child_seqs[ blockIdx.x * child_seuqnece_width + tid ];
    }
    __syncthreads();

    c = (( p & ~c ) | ( q & c ));

    if( tid < child_sequence_width ) {
        child_seqs[ blockIdx.x * child_sequence_width + tid ] = c;
    }
}

#endif  // RECOMBINER_CUH_
