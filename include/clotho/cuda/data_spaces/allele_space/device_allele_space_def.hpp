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

template < class RealType, class IntType, class OrderTag >
struct device_allele_space {
    typedef RealType    real_type;
    typedef IntType     int_type;

    static const unsigned int ALIGNMENT_SIZE = 128;

    real_type   * locations;
    int_type    * free_list;

    unsigned int    size, capacity;
    unsigned int    free_count;

    __host__ __device__ unsigned int allele_count( ) {
        return size;
    }

    __host__ __device__ unsigned int free_list_size( ) {
        return compute_free_list_size< IntType >( size );
    }
};


template < class RealType, class IntType, class OrderTag >
std::ostream & operator<<( std::ostream & out, const device_allele_space< RealType, IntType, OrderTag > & rhs ) {
    out << "{";
    out << "\"locations\": " << std::hex << rhs.locations;
    out << ",\n\"free_list\": " << std::hex << rhs.free_list;
    out << std::dec;
    out << ",\n\"free_count\": " << rhs.free_count;
    out << ",\n\"size\": " << rhs.size;
    out << ",\n\"capacity\": " << rhs.capacity;
    out << "}";
    return out;
}

#endif  // DEVICE_ALLELE_SPACE_DEF_HPP_
