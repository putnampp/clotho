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

template < class AlleleSpaceType, class IntType >
__global__ void _merge_space( AlleleSpaceType * in_space
                            , device_free_space< IntType, unit_ordered_tag< IntType > > * in_free
                            , device_event_space< IntType, unit_ordered_tag< IntType > > * evts
                            , device_free_space< IntType, unit_ordered_tag< IntType > > * out_free
                            , AlleleSpaceType * out_space ) {
    typedef AlleleSpaceType  allele_space_type;
    typedef device_event_space< IntType, unit_ordered_tag< IntType > >  event_space_type;
    
    typedef typename event_space_type::int_type         int_type;
    typedef typename event_space_type::order_tag_type   order_tag_type;

    if( (blockIdx.y * gridDim.x + blockIdx.x) != 0 ) return;

    assert( order_tag_type::OBJECTS_PER_UNIT == 32 );
    assert( (blockDim.x * blockDim.y) == order_tag_type::OBJECTS_PER_UNIT );

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    int_type lane_free = in_free->bin_summary[tid];
    int_type lane_count = evts->bin_summary[tid];

    int lane_overflow = lane_count - (( lane_free < lane_count ) ? lane_free : lane_count);    // determine whether there will be an overflow in this bin


    for( unsigned int i = 1; i < 32; i <<= 1 ) {    // move the max overflow to lane_id==31
        int t = __shfl_up( lane_overflow, i );
        lane_overflow = ((lane_overflow < t ) ? t : lane_overflow);
    }

    if( tid == 31 ) {
        // only one thread should call resize_space_impl
        // lane_overflow for thread 31 has the maximum overflow size for all bins
        //
        
        unsigned int size = in_space->capacity;
        unsigned int N = size + lane_overflow * order_tag_type::OBJECTS_PER_UNIT;
        _resize_space_impl( out_space, N );  // will pad space based upon ALIGNMENT_SIZE

        N = out_space->capacity;
        _resize_space_impl( out_free, N );

        if( size > 0 ) {
            _update_space( in_space, out_space );
            _update_space( in_free, out_free );
        }
    }
}

#endif  // DEVICE_ALLELE_SPACE_UNIT_ORDERED_KERNELS_HPP_
