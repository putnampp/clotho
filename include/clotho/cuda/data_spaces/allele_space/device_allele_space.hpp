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

#include "clotho/cuda/data_spaces/allele_space/device_allele_space_def.hpp"
#include "clotho/cuda/data_spaces/allele_space/device_allele_space_api.hpp"

#include "clotho/cuda/data_spaces/free_space/device_free_space.hpp"
#include "clotho/cuda/data_spaces/event_space/device_event_space.hpp"

#include "clotho/cuda/data_spaces/allele_space/device_allele_space_unordered_kernels.hpp"

#include "clotho/cuda/data_spaces/allele_space/device_allele_space_helper.hpp"
#include "clotho/cuda/data_spaces/allele_space/merge_space_helper.hpp"

template < class RealType/*, class IntType, class OrderTag*/ >
__device__ bool _resize_space_impl( device_allele_space< RealType/*, IntType, OrderTag*/ > * aspace, unsigned int N ) {
    typedef device_allele_space< RealType/*, IntType, OrderTag*/ > space_type;
    if( N % space_type::ALIGNMENT_SIZE != 0 ) {
        N -= (N % space_type::ALIGNMENT_SIZE);
        N += space_type::ALIGNMENT_SIZE;
    }

    bool space_increased = false;
    if( aspace->capacity < N ) {
        if( aspace->locations ) {
            delete aspace->locations;
        }
        aspace->locations = new typename space_type::real_type[N];
        
//        if( aspace->free_list ) {
//            delete aspace->free_list;
//        }

//        unsigned int free_size = compute_free_list_size< typename space_type::int_type >( N );
//        aspace->free_list = new typename space_type::int_type[free_size];

        memset(aspace->locations, 0, N * sizeof( typename space_type::real_type ) );
//        memset(aspace->free_list, -1, free_size * sizeof( typename space_type::int_type ) );

        aspace->capacity = N;
        space_increased = true;
    }

    aspace->size = N;
    return space_increased;
}

template < class RealType >
__device__ void _resize_space_impl( device_weighted_allele_space< RealType > * aspace, unsigned int N ) {
    typedef device_allele_space< RealType > space_type;

    if( _resize_space_impl( (device_allele_space< RealType > *) aspace, N ) ) {
        if( aspace->weights != NULL ) {
            delete aspace->weights;
        }

        aspace->weights = new typename space_type::real_type[ N ];
    }    
}

template < class RealType >
__global__ void _resize_space( device_allele_space< RealType > * aspace, unsigned int N ) {
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    if( tid == 0 ) {
        _resize_space_impl( aspace, N );
    }
}

template < class RealType >
__global__ void _resize_space( device_weighted_allele_space< RealType > * aspace, unsigned int N ) {
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    if( tid == 0 ) {
        _resize_space_impl( aspace, N );
    }
}

template < class RealType >
__global__ void _delete_space( device_allele_space< RealType > * aspace ) {
    typedef device_allele_space< RealType > space_type;

    space_type local = *aspace;

    if( local.locations ) {
        delete local.locations;
    }
}

template < class RealType >
__global__ void _delete_space( device_weighted_allele_space< RealType > * aspace ) {
    typedef device_weighted_allele_space< RealType > space_type;

    space_type local = *aspace;

    if( local.locations ) {
        delete local.locations;
        delete local.weights;
    }
}


template < class RealType >
void merge_allele_space( device_allele_space< RealType > * a
                        , device_allele_space< RealType > * b
                        , device_allele_space< RealType > * output ) {

}

template < class RealType, class IntType, class OrderTag >
void merge_space( device_allele_space< RealType/*, IntType, OrderTag*/ > * in_space
                , device_free_space< IntType, OrderTag > * fspace
                , device_event_space< IntType, OrderTag > * evts
                , device_allele_space< RealType/*, IntType, OrderTag*/ > * out_space ) {

    typedef merge_execution_config< OrderTag > config_type;
    _merge_space<<< config_type::BLOCK_COUNT, config_type::THREAD_COUNT >>>( in_space, fspace, evts, out_space );
}

#endif  // DEVICE_ALLELE_SPACE_HPP_
