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
#ifndef DEVICE_ALLELE_SPACE_HPP_
#define DEVICE_ALLELE_SPACE_HPP_

#include "clotho/cuda/allele_space/device_allele_space_def.hpp"
#include "clotho/cuda/allele_space/device_allele_space_api.hpp"

#include "clotho/cuda/allele_space/device_allele_space_unordered.hpp"
#include "clotho/cuda/allele_space/device_allele_space_partially_ordered.hpp"

template < class RealType, class IntType, class OrderTag >
void allele_space_alloc( device_allele_space< RealType, IntType, OrderTag > *& a, size_t N ) {
    size_t n = sizeof( device_allele_space< RealType, IntType, OrderTag > );
    assert( cudaMalloc( (void ** ) &a, n ) == cudaSuccess );
    assert( cudaMemset( a, 0, n ) == cudaSuccess );

    if( N ) {
        resize_space( a, N );
    }
}

template < class RealType, class IntType, class OrderTag >
void allele_space_free( device_allele_space< RealType, IntType, OrderTag > * a ) {
    typedef device_allele_space< RealType, IntType, OrderTag > space_type;

    space_type hCopy;
    assert( cudaMemcpy( &hCopy, a, sizeof(space_type), cudaMemcpyDeviceToHost ) == cudaSuccess );

    cudaFree( hCopy.locations );
    cudaFree( hCopy.free_list );

    cudaFree( a );
}

template < class RealType, class IntType, class OrderTag >
unsigned int resize_space( device_allele_space< RealType, IntType, OrderTag > * a, unsigned int N ) {
    typedef device_allele_space< RealType, IntType, OrderTag > space_type;
    space_type hSpace;
    cudaMemcpy( &hSpace, a, sizeof( space_type ), cudaMemcpyDeviceToHost );

    if( hSpace.capacity < N ) {
        // need to increase capacity
        std::cerr << "Increasing allele space capacity and size: " << hSpace.capacity << " -> " << N << std::endl;
        typename space_type::real_type * locs;

        size_t loc_size = N * sizeof( typename space_type::real_type );
        assert( cudaMalloc( (void **) &locs, loc_size) == cudaSuccess );

        typename space_type::int_type * free;
        size_t free_size = compute_free_list_size< typename space_type::int_type>( N ) * sizeof( typename space_type::int_type );
        assert( cudaMalloc( (void **) &free, free_size) == cudaSuccess );
        assert( cudaMemset( free, 255, free_size ) == cudaSuccess );

        if( hSpace.allele_count() ) {
            // Copy existing lists
            assert( cudaMemcpy( locs, hSpace.locations, hSpace.allele_count() * sizeof( typename space_type::real_type), cudaMemcpyDeviceToDevice ) == cudaSuccess );
            assert( cudaMemcpy( free, hSpace.free_list, hSpace.free_list_size() * sizeof( typename space_type::int_type), cudaMemcpyDeviceToDevice ) == cudaSuccess);

            cudaFree( hSpace.locations );
            cudaFree( hSpace.free_list );
        }

        hSpace.locations = locs;
        hSpace.free_list = free;
        hSpace.size = N;
        hSpace.capacity = N;

        cudaMemcpy( a, &hSpace, sizeof( space_type ), cudaMemcpyHostToDevice );

        update_free_count( a );
    } else if( hSpace.size != N ) {
        std::cerr << "Resizing allele space: " << hSpace.size << " -> " << N << std::endl;
        hSpace.size = N;

        cudaMemcpy( a, &hSpace, sizeof( space_type ), cudaMemcpyHostToDevice );
        update_free_count( a );
    }

    return N;
}

template < class RealType, class IntType, class OrderTag >
void update_free_count( device_allele_space< RealType, IntType, OrderTag > * a ) {
    _update_free_count<<< 1, 32 >>>( a );
}

#endif  // DEVICE_ALLELE_SPACE_HPP_
