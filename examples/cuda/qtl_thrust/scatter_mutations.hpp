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
#ifndef SCATTER_MUTATIONS_HPP_
#define SCATTER_MUTATIONS_HPP_

#include <cuda.h>

__global__ void scatter_mutations( double * rands
                                 , double * alleles
                                 , unsigned int * free_list
                                 , unsigned int * sequences
                                 , unsigned int * mut_events
                                 , unsigned int seq_count
                                 , unsigned int allele_count
                                );

__global__ void scatter_mutations( double * rands
                                 , double * alleles
                                 , unsigned int * free_list
                                 , unsigned int * sequences
                                 , unsigned int * mut_events
                                 , unsigned int seq_count
                                 , unsigned int allele_count
                                 , unsigned int * dbg
                                );

#endif  // SCATTER_MUTATIONS_HPP_
