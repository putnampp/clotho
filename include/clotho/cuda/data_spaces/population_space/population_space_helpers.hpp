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
#ifndef POPULATION_SPACE_HELPERS_HPP_
#define POPULATION_SPACE_HELPERS_HPP_

#include "clotho/cuda/data_spaces/population_space/population_space_helper_kernels.hpp"

template < class IntType, class OrderTag >
void update_free_map( device_free_space< IntType, OrderTag > * free_space ) {
    dim3 blocks(1,1,1), threads( 32, 1, 1 );
    update_free_map_kernel<<< blocks, threads >>>( free_space );
}

template < class IntType, class OrderTag >
void validate_free_map( device_free_space< IntType, OrderTag > * free_space, device_event_space< IntType, OrderTag > * events ) {
    validate_free_space_size<<< 1, 1 >>>( free_space, events );
}

/*
template < class IntType, class OrderTag >
void update_free_space( device_sequence_space< IntType > * seq_space
                       , device_free_space< IntType, OrderTag > * free_space ) {
    resize_free_space_kernel<<< 1, 1 >>>( seq_space, free_space);

    update_free_space_kernel2<<< 200, 32 >>>( seq_space, free_space );

    clear_free_space_kernel2<<< 200, 32 >>> ( seq_space, free_space );

    update_free_map( free_space );

    cudaDeviceSynchronize();
}*/
template < class IntType, class OrderTag >
void update_free_space2( device_sequence_space< IntType > * seq_space
                            , device_free_space< IntType, OrderTag > * free_space ) {
    clotho::utility::algo_version< 1 > * v2 = NULL;
    update_free_space_kernel<<< 200, 32 >>>( seq_space, free_space, v2 );

    dim3 blocks( 2, 1, 1), threads( 32, 1, 1);
    update_free_space_total_kernel<<< blocks, threads >>>( free_space );

    clotho::utility::algo_version< 4 > * v = NULL;
    clear_free_space_kernel<<< 10, 1024 >>> ( seq_space, free_space, v );
}

template < class RealType, class IntType, class OrderTag >
std::ostream & operator<<( std::ostream & out, const PopulationSpace< RealType, IntType, OrderTag > & pop ) {
    out << "{";
    out << "\n\"sequences\": " << pop.sequences;
    out << ",\n\"alleles\": " << pop.alleles;
    out << ",\n\"free_space\": " << pop.free_space;
    out << "\n}";

    return out;
}
#endif  // POPULATION_SPACE_HELPERS_HPP_
