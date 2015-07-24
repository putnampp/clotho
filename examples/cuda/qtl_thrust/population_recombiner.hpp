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
#ifndef POPULATION_RECOMBINER_HPP_
#define POPULATION_RECOMBINER_HPP_

#include <cuda.h>

extern const unsigned int SEQ_PER_INDIVIDUAL;
extern const unsigned int PARENT_PER_OFF;
extern const unsigned int RANDOM_PER_PARENT;   // id offset & swapping of source sequences
extern const unsigned int MAX_OFFSPRING;
extern const unsigned int MAX_BLOCKS_PER_STRIDE;

extern const unsigned int MAX_OFFSPRING_SEQ;    // 16
extern const unsigned int MAX_PARENTS;          // 16
extern const unsigned int MAX_PARENT_SEQ;       // 32

extern const unsigned int MAX_PARENT_BLOCKS;    // 32 * 32 = 1024
extern const unsigned int MAX_OFFSPRING_BLOCKS; // 16 * 32 = 512
extern const unsigned int MAX_RANDOMS;          // 16 * 2 = 32

extern const unsigned int MAX_THREADS;          // 512

__global__ void recombine_population( double * rand_pool
                                    , unsigned int * parents
                                    , unsigned int * offspring
                                    , unsigned int parent_rows
                                    , unsigned int parent_cols
                                    , unsigned int off_rows
                                    , unsigned int off_cols);


__global__ void recombine( unsigned int * p0
                            , unsigned int * p1
                            , unsigned int * off
                            , unsigned int cols );

__global__ void recombiner( double * rands
                            , unsigned int * parents
                            , unsigned int parent_rows
                            , unsigned int parent_cols
                            , unsigned int * off
                            , unsigned int cols
                            , unsigned int seq_offset );

#endif  // POPULATION_RECOMBINER_HPP_
