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
#ifndef DEVICE_ALLELE_SPACE_DEF_HPP_
#define DEVICE_ALLELE_SPACE_DEF_HPP_

#include <cuda.h>

template < class IntType >
static __host__ __device__ unsigned int compute_free_list_size( unsigned int N ) {
    return (N / (sizeof(IntType) * 8));
}

template < class RealType >
struct device_allele_space {
    typedef RealType    real_type;

    static const unsigned int ALIGNMENT_SIZE = 128;

    real_type       * locations;
    real_type       * neutral;

    real_type       neutral_p;
    unsigned int    size, capacity;
};

template < class RealType >
struct device_weighted_allele_space : public device_allele_space< RealType > {
    typename device_allele_space< RealType >::real_type   * weights;
};

#endif  // DEVICE_ALLELE_SPACE_DEF_HPP_
