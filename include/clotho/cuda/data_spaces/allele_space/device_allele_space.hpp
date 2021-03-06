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
#include "clotho/cuda/data_spaces/allele_space/device_allele_space_unit_ordered_kernels.hpp"

#include "clotho/cuda/data_spaces/allele_space/device_allele_space_helper.hpp"
#include "clotho/cuda/data_spaces/allele_space/merge_space_helper.hpp"
#include "clotho/cuda/helper_macros.hpp"

#include "clotho/cuda/popcount_kernel.h"

template < class RealType >
__device__ bool _resize_space_impl( device_allele_space< RealType > * aspace, unsigned int N, bool copy_content = false ) {
//    assert( blockIdx.y * gridDim.x + blockIdx.x == 0 );
//    assert( threadIdx.y * blockDim.x + threadIdx.x == 0);

    typedef device_allele_space< RealType > space_type;
    typedef typename space_type::real_type  real_type;

    unsigned int tail = N % space_type::ALIGNMENT_SIZE;
    if( tail != 0 ) {
        N -= tail;
        N += space_type::ALIGNMENT_SIZE;
    }

    unsigned int M = aspace->capacity;
    bool cap_increased = false;
    if( M < N ) {
        real_type * l = new real_type[ N ];
        real_type * neu = new real_type[ N ];

        assert( l != NULL );

        memset(l, 0, N * sizeof( real_type ) );
        if( aspace->locations ) {
            if( copy_content ) {
                memcpy(l, aspace->locations, M * sizeof(real_type));
                memcpy(neu, aspace->neutral, M * sizeof(real_type));
            }
            delete aspace->locations;
            delete aspace->neutral;
        }


        aspace->locations = l;
        aspace->neutral = neu;
        aspace->capacity = N;
        cap_increased = true;
    }

    aspace->size = N;

    return cap_increased;
}

template < class RealType >
__device__ bool  _resize_space_impl( device_weighted_allele_space< RealType > * aspace, unsigned int N, bool copy_content = false ) {
    typedef device_allele_space< RealType > space_type;
    typedef typename space_type::real_type  real_type;

    unsigned int M = aspace->capacity;
    bool cap_increase =  _resize_space_impl( (device_allele_space< RealType > *) aspace, N, copy_content );

    if( cap_increase ) {
        N = aspace->capacity;

        real_type * w = new typename space_type::real_type[ N ];
        assert( w != NULL );

        memset( w, 0, N * sizeof( real_type ) );
        if( aspace->weights != NULL ) {
            if( copy_content ) {
                memcpy(w, aspace->weights, M * sizeof( real_type ));    
            }
        
            delete aspace->weights;
        }

        aspace->weights = w;
    }

    return cap_increase;
}

template < class RealType >
__global__ void _resize_space( device_allele_space< RealType > * aspace, unsigned int N ) {
    assert( blockIdx.y * gridDim.x + blockIdx.x == 0 );

    if( threadIdx.y * blockDim.x + threadIdx.x == 0 ) {
        _resize_space_impl( aspace, N );
    }
}

template < class RealType >
__global__ void _resize_space( device_weighted_allele_space< RealType > * aspace, unsigned int N ) {
    assert( blockIdx.y * gridDim.x + blockIdx.x == 0 );

    if( threadIdx.y * blockDim.x + threadIdx.x == 0 ) {
        _resize_space_impl( aspace, N );
    }
}

template < class RealType >
__device__ void _delete_space_impl( device_allele_space< RealType > * aspace ) {
    if( blockIdx.y * gridDim.x + blockIdx.x > 0 ) return;

    if( aspace->locations  && threadIdx.y * blockDim.x + threadIdx.x == 0) {
        delete aspace->locations;
        aspace->locations = NULL;

        delete aspace->neutral;
        aspace->neutral = NULL;
    }
}

template < class RealType >
__device__ void _delete_space_impl( device_weighted_allele_space< RealType > * aspace ) {
    if( blockIdx.y * gridDim.x + blockIdx.x > 0 ) return;

    _delete_space_impl( (device_allele_space< RealType > * ) aspace);

    if( aspace->weights && threadIdx.y * blockDim.x + threadIdx.x == 0 ) {
        delete aspace->weights;
        aspace->weights = NULL;
    }
}

template < class AlleleSpaceType >
__global__ void _delete_space( AlleleSpaceType * aspace ) {
    
    if( blockIdx.y * gridDim.x + blockIdx.x > 0 ) return;

    if( threadIdx.y * blockDim.x + threadIdx.x == 0 ) {
        _delete_space_impl( aspace );
    }
}

template < class RealType >
void merge_allele_space( device_allele_space< RealType > * a
                        , device_allele_space< RealType > * b
                        , device_allele_space< RealType > * output ) {

}

template < class AlleleSpaceType, class IntType, class OrderTag >
void merge_space( AlleleSpaceType * in_space
                , device_free_space< IntType, OrderTag > * fspace
                , device_event_space< IntType, OrderTag > * evts
                , device_free_space< IntType, OrderTag > * ofspace
                , AlleleSpaceType * out_space ) {

    typedef merge_execution_config< OrderTag > config_type;
    expand_spaces_kernel<<< config_type::BLOCK_COUNT, config_type::THREAD_COUNT >>>( in_space, fspace, evts, ofspace, out_space );
    CHECK_LAST_KERNEL_EXEC

    update_space_kernel<<< 3, 1024 >>>( in_space, out_space );
    CHECK_LAST_KERNEL_EXEC

    update_space_kernel<<< 3, 1024 >>>( fspace, ofspace );
    CHECK_LAST_KERNEL_EXEC
}

/*
template < class RealType, class IntType, class OrderTag >
void merge_space( device_weighted_allele_space< RealType > * in_space
                , device_free_space< IntType, OrderTag > * fspace
                , device_event_space< IntType, OrderTag > * evts
                , device_free_space< IntType, OrderTag > * ofspace
                , device_weighted_allele_space< RealType > * out_space ) {

    typedef merge_execution_config< OrderTag > config_type;

    _merge_space<<< config_type::BLOCK_COUNT, config_type::THREAD_COUNT >>>( in_space, fspace, evts, ofspace, out_space );
    CHECK_LAST_KERNEL_EXEC
}*/

template < class RealType >
__device__ bool _move_allele( device_allele_space< RealType > * in_space
                                , unsigned int in_idx
                                , device_allele_space< RealType > * out_space
                                , unsigned int out_idx ) {
    unsigned int in_size = in_space->capacity;
    unsigned int out_size = out_space->capacity;

    bool success = ((in_idx < in_size) && (out_idx < out_size));

    if( success ) {
        out_space->locations[out_idx] = in_space->locations[in_idx];
        out_space->neutral[out_idx] = in_space->neutral[in_idx];
    }   

    return success;
}

template < class RealType >
__device__ bool _move_allele( device_weighted_allele_space< RealType > * in_space
                                , unsigned int in_idx
                                , device_weighted_allele_space< RealType > * out_space
                                , unsigned int out_idx ) {
    typedef device_allele_space< RealType > base_type;

    bool success = _move_allele( (base_type *) in_space, in_idx, (base_type *) out_space, out_idx);
    if( success ) {
        out_space->weights[ out_idx ] = in_space->weights[ in_idx ];
    }
    return success;
}

template < class AlleleSpaceType, class IntType, class OrderTag >
__global__ void resize_fixed_allele_kernel( AlleleSpaceType * alls, device_free_space< IntType, OrderTag > * free_space ) {

    assert( blockIdx.y * gridDim.x + blockIdx.x == 0 );

    unsigned int _count = free_space->fixed_count;
    if( _count == 0 ) return;

    //printf( "%d fixed allele encountered\n", _count );
    //


    if( threadIdx.y * blockDim.x + threadIdx.x == 0 ) {
        _count += alls->size;

        _resize_space_impl( alls, _count, true );
    }
}

template < class AlleleSpaceType, class IntType, class OrderTag >
__global__ void move_fixed_allele_kernel( AlleleSpaceType * dest, AlleleSpaceType * src, device_free_space< IntType, OrderTag > * src_free_space ) {
    unsigned int _count = src_free_space->fixed_count;
    if( _count == 0 ) {
        return;
    }

    assert( blockIdx.y * gridDim.x + blockIdx.y == 0 );   // algorithm assumes single block

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int lane_id = (tid & 31);

    unsigned int N = dest->size;

    if( tid == 0 ) {
        _resize_space_impl( dest, N + _count, true );
    }
    __syncthreads();

    typedef device_free_space< IntType, OrderTag >  space_type;
    typedef typename space_type::int_type           int_type;

    assert( sizeof(int_type) * 8 == space_type::OBJECTS_PER_INT );

    int_type * xlist = src_free_space->fixed_list;
    unsigned int M = src_free_space->size;

    unsigned int offset = 0, a_idx = lane_id;
    unsigned int lane_mask = (1 << lane_id);
    unsigned int lane_lower_mask = (lane_mask - 1);

    popcountGPU< int_type > pc;

    _count += N;

    while( N < _count && offset < M ) {
        int_type x = xlist[offset++];

        if( x != 0 ) {  // true for all threads in block
            unsigned int offset = pc.evalGPU( (x & lane_lower_mask) );

            if( (x & lane_mask) ) {
                N += offset;
                _move_allele(src, a_idx, dest, N );
                ++N;
            }
            __syncthreads();

            for( unsigned int i = 1; i < 32; i <<= 1 ) {
                unsigned int n = __shfl_up( N, i );
                N = ((n > N) ? n : N);
            }

            N = __shfl( N, 31 );
        }

        a_idx += space_type::OBJECTS_PER_INT;
    }

    if( lane_id == 0 ) {
        dest->size = N;
    }

//    assert( N >= _count);
}

template < class StateType, class RealType >
__device__ RealType _generate_random_allele_location( StateType * state, RealType p = 0.0) {
    RealType x = 0.0;
    do {
        x = curand_uniform( state );   // curand_uniform in (0,1]
    } while( x >= 1.0);

    return x;
}

template < class StateType, class RealType >
__device__ RealType _generate_neutrality( StateType * state, RealType p = 0.0 ) {
    if( p == 0.0) return 0.0;
    if( p == 1.0) return 1.0;

    RealType x = curand_uniform( state );

    // if uniform random number x is less than p
    // then allele is neutral represented by scaling
    // factor of 0.0
    x = ((x <= p) ? 0.0 : 1.0);

    return x;
}

template < class StateType, class RealType >
__device__ void _generate_random_allele( StateType * state
                                        , device_allele_space< RealType > * alleles
                                        , unsigned int idx 
                                        , unordered_tag * tag ) {
    RealType x = _generate_random_allele_location( state, (RealType) 0.0 );
    alleles->locations[ idx ] = x;  // x in (0,1) == mutation cannot occur at very beginning or end of a sequence

    RealType p = alleles->neutral_p;

    RealType n = _generate_neutrality( state, p );

    alleles->neutral[ idx ] = n;
}

template < class StateType, class RealType, class IntType >
__device__ void _generate_random_allele( StateType * state
                                        , device_allele_space< RealType > * alleles
                                        , unsigned int idx 
                                        , unit_ordered_tag< IntType > * tag ) {
    
    RealType x = _generate_random_allele_location( state );

    RealType lane_id = (RealType )(idx & 31);

    x = (( lane_id + x ) / ((RealType) unit_ordered_tag< IntType >::OBJECTS_PER_UNIT));

    RealType p = alleles->neutral_p;
    RealType n = _generate_neutrality< RealType >( state, p );

    alleles->locations[ idx ] = x;
    alleles->neutral[ idx ] = n;
}

template < class StateType, class RealType, class OrderTag >
__device__ void _generate_random_allele( StateType * state
                                        , device_weighted_allele_space< RealType > * alleles
                                        , unsigned int idx 
                                        , OrderTag * tag ) {
    _generate_random_allele( state, (device_allele_space< RealType > * ) alleles, idx, tag );

    RealType n = alleles->neutral[ idx ];

    RealType y = n * curand_normal( state );    // 0 weights for neutral alleles; simplifies phenotype calculation algorithm
    alleles->weights[ idx ] = y;
}

#endif  // DEVICE_ALLELE_SPACE_HPP_
