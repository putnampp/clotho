#include "square.h"

#include <stdio.h>
#include <cuda.h>

__global__ void square( unsigned int * a, int N ) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if( idx < N ) a[idx] = idx*idx;
}

Square::Square() {}
Square::~Square() {}

void Square::operator()() {
    unsigned int * dest;

    size_t size = N * sizeof(unsigned int);
    cudaMalloc( (void **) &dest, size);

    square<<<1,N>>>( dest, N );

    cudaMemcpy( m_a, dest, size, cudaMemcpyDeviceToHost );
    cudaFree( dest );
}

std::ostream & operator<<( std::ostream & out, const Square & rhs ) {
    for( unsigned int i = 0; i < Square::N; ++i ) {
        out << i << " -> " << rhs.m_a[i] << std::endl;
    }
    return out;
}
