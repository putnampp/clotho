#include "square.h"

#include <stdio.h>
#include <cuda.h>
#include <iostream>
#include <curand_kernel.h>

__global__ void initRNG( curandState * rngStates, const unsigned int seed ) {
    unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

    curand_init( seed, idx, 0, &rngStates[idx] );
}

__global__ void square( Square::int_type * a, int N ) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if( idx < N ) a[idx] = idx * idx;
}

__global__ void squareRNG( Square::int_type * a, int N, curandState * rngStates ) {

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    curandState lState = rngStates[idx];

    Square::int_type r = curand(&lState);
    if( idx < N ) {
        a[idx] = r;
    }
}

Square::Square() : m_a(NULL), m_size(1), m_maxBlocks(1), m_maxThreadsPerBlock(256) {
    init();
}

void Square::init() {
    cudaDeviceProp m_cdp;
    cudaError_t err = cudaGetDeviceProperties( &m_cdp, 0 );

    if( err != cudaSuccess ) {
        std::cerr << "Unable to get device properties" << std::endl;
    } else {
        std::cerr << "Maximum Threads Per Block: " << m_cdp.maxThreadsPerBlock << std::endl;
        m_maxThreadsPerBlock = m_cdp.maxThreadsDim[0];
        m_maxBlocks = m_cdp.maxGridSize[0];
    }

//    std::cerr << "Maximum Threads Per Block: " << m_maxThreadsPerBlock << std::endl;
//    std::cerr << "Maximum Blocks: " << m_maxBlocks << std::endl;

    m_size = ((m_maxBlocks > 6 * m_maxThreadsPerBlock) ? 6 * m_maxThreadsPerBlock * m_maxThreadsPerBlock : m_maxBlocks * m_maxThreadsPerBlock);

    m_a = new int_type[ m_size ];

    size_t size = m_size * sizeof(int_type);
    cudaMalloc( (void **) &m_dest, size);

//    std::cerr << "Compute Mode: " <<  m_cdp.computeMode << std::endl;
//    std::cerr << "Device Overlap: " <<  m_cdp.deviceOverlap << std::endl;
}

Square::~Square() {
    if( m_a ) delete [] m_a;
    cudaFree( m_dest );
}

size_t Square::size() const { return m_size; }

void Square::operator()() {
    int block_count = (m_size / m_maxThreadsPerBlock );
    square<<<block_count,m_maxThreadsPerBlock>>>( m_dest, m_size );

    cudaMemcpy( m_a, m_dest, m_size * sizeof(int_type), cudaMemcpyDeviceToHost );
}

void Square::random_list() {
    int block_count = (m_size / m_maxThreadsPerBlock );
    curandState *d_rngStates = 0;
    cudaMalloc( (void **) &d_rngStates, block_count * sizeof( d_rngStates ) );

    initRNG<<< 1, block_count >>>( d_rngStates, 1234 );

    squareRNG<<<1, block_count,m_maxThreadsPerBlock>>>( m_dest, m_size, d_rngStates );

    cudaMemcpy( m_a, m_dest, m_size * sizeof(int_type), cudaMemcpyDeviceToHost );
    cudaFree( d_rngStates );
}

std::ostream & operator<<( std::ostream & out, const Square & rhs ) {
    for( unsigned int i = 0; i < rhs.size(); ++i ) {
        out << i << " -> " << rhs.m_a[i] << std::endl;
    }
    return out;
}
