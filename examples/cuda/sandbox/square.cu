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
#include "square.h"

#include <stdio.h>
#include <iostream>
#include <cassert>

// from CUDA documentation
//const unsigned int maxThreadsPerState = 256;

// each thread initializes a random state
__global__ void initRNG( curandState * rngStates, Square::seed_type * seeds, unsigned int N ) {
    unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    unsigned int tmax = blockDim.x * blockDim.y;

    unsigned int idx = bid * tmax + tid;

    __shared__ unsigned int tmp[ 1024 ];
    if( idx < N ) {
        tmp[ tid ] = seeds[ idx ];
        __syncthreads();

        curandState localState;
        curand_init( tmp[tid], tid, 0, &localState  );
        rngStates[ idx ] = localState;
    }
}

/*
__global__ void initRNG( curandStateMtgp32_t * rngStates, Square::seed_type * seeds, unsigned int N ) {
    unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

    curand_init( seed, idx, 0, &rngStates[idx] );
}*/

__global__ void square( Square::int_type * a, int N ) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if( idx < N ) a[idx] = idx * idx;
}

template < class State >
__global__ void squareRNG( Square::int_type * a, int N, State * rngStates ) {

    unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    unsigned int idx = bid * (blockDim.x * blockDim.y) + tid;
    __shared__ Square::int_type tmp[ 1024 ];    // assume max 1024 threads/block

    // every thread generates a random number
    if( idx < N ) {
        tmp[tid] = curand(rngStates + (bid * blockDim.y  + threadIdx.y) );
        a[ idx ] = tmp[tid] * tmp[tid];    // copy square to global memory
    }
}

Square::Square( boost::random::mt19937 & rng ) :
    m_a(NULL)
    , m_dest(NULL)
    , m_size(0)
    , m_capacity(0)
    , m_maxBlocks(0)
    , m_maxThreadsPerBlock(0)
    , m_status(true)
    , m_dStates( &rng )
{
    init();
}

void Square::init() {
    cudaDeviceProp m_cdp;
    cudaError_t err = cudaGetDeviceProperties( &m_cdp, 0 );

    if( err != cudaSuccess ) {
        std::cerr << "Unable to get device properties" << std::endl;
        m_status = false;
        return;
    } else {
        std::cerr << "Maximum Threads Per Block: " << m_cdp.maxThreadsPerBlock << std::endl;
        m_maxThreadsPerBlock = m_cdp.maxThreadsPerBlock;
        m_maxBlocks = m_cdp.maxGridSize[0];
    }


    // initialize random states
//    unsigned int state_count = m_maxThreadsPerBlock / maxThreadsPerState;
//
//    if( m_maxThreadsPerBlock % maxThreadsPerState ) { ++state_count; }
//    m_dStates.resize( state_count );
}

Square::~Square() {
    if( m_a ) free(m_a);
    if( m_dest ) cudaFree( m_dest );
}

size_t Square::size() const { return m_size; }

bool Square::good() const {
    return m_status && m_dStates.good();
}

void Square::operator()( unsigned int s ) {
    if( !good() ) return;

    resize( s );

    assert( s == m_size );

    int block_count = (m_size / m_maxThreadsPerBlock );
    if( m_size % m_maxThreadsPerBlock ) { ++block_count; }

    dim3 gdim( block_count, 1, 1), 
         bdim = curand_state_type::makeBlockDimension( m_maxThreadsPerBlock );

    m_dStates.updateStates( gdim.x * gdim.y * bdim.y );

//    std::cerr << "Block Dimensions: <" << bdim.x << ", " << bdim.y << ", " << bdim.z << ">" << std::endl;
//    std::cerr << "Grid Dimensions: <" << gdim.x << ", " << gdim.y << ", " << gdim.z << ">" << std::endl;

    squareRNG<<< gdim, bdim >>>( m_dest, m_size, m_dStates.getStates() );

    //int_type * d = m_dest;
    //curand_state_type::pointer pS = m_dStates.getStates();
    //while( s ) {
    //    if( s >= bdim.x ) {
    //        squareRNG<<< 1, bdim.x >>>( d, bdim.x, pS );
    //        d += bdim.x;
    //        s -= bdim.x;
    //        ++pS;
    //    } else {
    //        squareRNG<<< 1, bdim.x >>>( d, s, pS );
    //        s = 0;
    //    }
    //}
    cudaMemcpy( m_a, m_dest, m_size * sizeof(int_type), cudaMemcpyDeviceToHost );
}

void Square::resize( unsigned int s ) {
    if( !good() ) return;

    if ( s > m_capacity ) {
        if( m_a ) {
            free( m_a );
        }
        if( m_dest ) {
            cudaFree( m_dest );
        }

        size_t byte_size = s * sizeof(int_type );
        m_a = (int_type *) malloc( byte_size );
        cudaError_t err = cudaMalloc( (void **) &m_dest, byte_size );
        if( err != cudaSuccess ) {
            std::cerr << "Unable to allocate device memory: " << std::endl;
        }

        m_capacity = s;
    }

    m_size = s;
}

/*
void Square::random_list() {
    int block_count = (m_size / m_maxThreadsPerBlock);

    if( m_size % m_maxThreadsPerBlock ) {
        ++block_count;
    }

    // max 256 threads/curandState
    unsigned int state_count = block_count * ( m_maxThreadsPerBlock / 256 );
    unsigned int seed_size = state_count * sizeof(unsigned int);
    unsigned int state_size = state_count * sizeof( curandState );

    unsigned int * seeds;
    seeds = malloc( seed_size );

    unsigned int * d_seeds;
    cudaMalloc( (void **) &d_seeds, seed_size );

    unsigned int t = state_count;
    while( t-- ) {
        seeds[ t ] = m_dist( *m_rng );
    }
 
    cudaMemcpy( d_seeds, seeds, seed_size, cudaMemcpyHostToDevice );
   
    curandState *d_rngStates = 0;
    cudaMalloc( (void **) &d_rngStates, state_size );

    unsigned int bcount = state_count / m_maxThreadsPerBlock;
    if( state_count % m_maxThreadsPerBlock ) { ++bcount; }

    unsigned int bx = bcount, by = 1;

    if( bcount > m_maxBlocks ) {
        by = bcount / m_maxBlocks;

        if( bcount % m_maxBlocks ) { ++by; }

        bx = m_maxBlocks;
    }

    unsigned int tx = m_maxThreadsPerBlock, ty = 1;


    initRNG<<< dim3(bx, by, 1), dim3(tx, ty, 1) >>>( d_rngStates, d_seeds, state_count );

    squareRNG<<<1, block_count, m_maxThreadsPerBlock >>>( m_dest, m_size, d_rngStates );

    free(seeds);    // should be performed on host during squareRNG on device

    cudaMemcpy( m_a, m_dest, m_size * sizeof(int_type), cudaMemcpyDeviceToHost );
    cudaFree( d_rngStates );
    cudaFree( d_seeds );
}*/

std::ostream & operator<<( std::ostream & out, const Square & rhs ) {
    if( rhs.good() ) {
        for( unsigned int i = 0; i < rhs.size(); ++i ) {
            out << i << " -> " << rhs.m_a[i] << "\n";
        }
    } else {
        out << "BAD STATE" << "\n";
    }
    return out;
}
