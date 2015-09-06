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
#ifndef DEVICE_ALLELE_SPACE_UNIT_ORDERED_HPP_
#define DEVICE_ALLELE_SPACE_UNIT_ORDERED_HPP_

#include "clotho/cuda/data_spaces/allele_space/device_allele_space_def.hpp"

#include "clotho/cuda/data_spaces/tags/unit_ordered_tag.hpp"

template < class RealType, class IntType >
struct device_allele_space< RealType, IntType, unit_ordered_tag< IntType > > {
    typedef RealType                    real_type;
    typedef IntType                     int_type;
    typedef unit_ordered_tag< IntType > order_tag_type;

    typedef unit_ordered_tag< IntType > order_type;

    static const unsigned int ALIGNMENT_SIZE = 128;

    real_type   * locations;
    int_type    * free_list;

    unsigned int free_count[ order_type::OBJECTS_PER_UNIT ];
    unsigned int size, capacity;

    __host__ __device__ unsigned int allele_count( ) {
        return size;
    }

    __host__ __device__ unsigned int free_list_size( ) {
        return compute_free_list_size< IntType >( size );
    }
};

/*
template < class RealType, class IntType >
std::ostream & operator<<( std::ostream & out, const device_allele_space< RealType, IntType, unit_ordered_tag< IntType >  > & rhs ) {
    out << "{";
    out << "\"locations\": " << std::hex << rhs.locations;
    out << ",\n\"free_list\": " << std::hex << rhs.free_list;

    out << std::dec;

    out << ",\n\"free_count\": [" << rhs.free_count[0];
    for( unsigned int i = 1; i < 32; ++i ) {
        out << ",\n  " << rhs.free_count[i];
    }

    out << "],\n\"size\": " << rhs.size;
    out << ",\n\"capacity\": " << rhs.capacity;
    out << "}";
    return out;
}*/

#endif  // DEVICE_ALLELE_SPACE_UNIT_ORDERED_HPP_
