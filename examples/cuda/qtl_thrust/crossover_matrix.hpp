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
#ifndef CROSSOVER_MATRIX_HPP_
#define CROSSOVER_MATRIX_HPP_

#include <cuda.h>

/**
 *  kernel expected execution dimensions  <<< dim3(pop_size / ( gpu blocking ability ), allele_count/GPU_THREADS (+ 1), 1), dim3(GPU_THREADS/32,32,1) >>>
 *
 *  max_dims.x = maximum number of alleles
 *  max_dims.y = maximum number of sequences
 *  max_dims.z = maximum number of random numbers
 */
__global__ void generate_crossover_matrix( double * rand_pool
                                         , double * alleles
                                         , unsigned int * event_list
                                         , unsigned int * sequences
                                         , dim3 max_dims );

__global__ void init_alleles( double * alleles, unsigned int allele_count );
__global__ void crossover2( unsigned int * seq, double rpoint );

__global__ void crossover( unsigned int * seq
                            , double * alleles
                            , unsigned int allele_count
                            , double rpoint);

#endif  // CROSSOVER_MATRIX_HPP_
