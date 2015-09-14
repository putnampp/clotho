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
#ifndef DEVICE_ALLELE_SPACE_UNIT_ORDERED_KERNELS_HPP_
#define DEVICE_ALLELE_SPACE_UNIT_ORDERED_KERNELS_HPP_

#include "clotho/cuda/data_spaces/allele_space/device_allele_space_unit_ordered.hpp"
#include "clotho/cuda/data_spaces/allele_space/device_allele_space_kernel_api.hpp"

/*
template < class RealType, class IntType >
__global__ void _update_free_count( device_allele_space< RealType, IntType, unit_ordered_tag< IntType > > * aspace ) {
    IntType tid = threadIdx.y *blockDim.x + threadIdx.x;
    IntType lane_mask = (1 << (tid & 31));

    IntType * local_free_list = aspace->free_list;
    unsigned int * fcounts = aspace->free_count;

    unsigned int N = aspace->free_list_size();
    unsigned int _count = 0;

    while( N-- ) {
        IntType f = *local_free_list;

        _count += !!(f & lane_mask);
        ++local_free_list;
    }

    fcounts[tid] = _count;
}*/

template < class RealType, class IntType >
__global__ void _merge_space( device_allele_space< RealType/*, IntType, unit_ordered_tag< IntType >*/ > * in_space
                            , device_free_space< IntType, unit_ordered_tag< IntType > > * fspace
                            , device_event_space< IntType, unit_ordered_tag< IntType > > * evts
                            , device_allele_space< RealType/*, IntType, unit_ordered_tag< IntType >*/ > * out_space ) {
    typedef device_allele_space< RealType/*, IntType, unit_ordered_tag< IntType >*/ >   allele_space_type;
    typedef device_event_space< IntType, unit_ordered_tag< IntType > >              event_space_type;
    
    typedef typename allele_space_type::int_type    int_type;
    typedef typename allele_space_type::real_type   real_type;

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    int_type in_free = fspace->free_count[tid];
    int_type bin = evts->bin_summary[tid];

    int n_units = 0;
    if( in_free < bin ) {   // determine whether there will be an overflow in this bin
        int tmp = bin - in_free;
        if( tmp > n_units ) {
            n_units = tmp;
        }
    }
    __syncthreads();

    for( unsigned int i = 1; i < 32; i <<= 1 ) {
        int t = __shfl_up( n_units, i );
        n_units = (( (tid >= i ) && (n_units < t) ) ? t : n_units);
    }

    if( tid == 31 ) {
        // only one thread should call resize_space_impl
        // n_units for thread 31 has the maximum overflow size for all bins
        //
        
        int_type size = in_space->size;
        int_type N = size + n_units * unit_ordered_tag< IntType >::OBJECTS_PER_UNIT;
        _resize_space_impl( out_space, N );  // will pad space based upon ALIGNMENT_SIZE

        real_type * locs = in_space->locations;
        
        if( locs ) {
            memcpy( out_space->locations, locs, size * sizeof( real_type ) );

//            size = compute_free_list_size< int_type >( size );
//            memcpy( out_space->free_list, in_space->free_list, size * sizeof( int_type ));
        }
    }
}

#endif  // DEVICE_ALLELE_SPACE_UNIT_ORDERED_KERNELS_HPP_
