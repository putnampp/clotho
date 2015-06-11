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
#include "square_host.h"
#include <stdint.h>
#include "popcount_kernel.h"

#include <iostream>

__global__ void square_mem( SquareHost::int_type * a, unsigned int N ) {
    unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    unsigned int idx = bid * (blockDim.x * blockDim.y) + tid;
    if( idx < N ) {
        a[idx] *= a[idx];
    }
}

SquareHost::SquareHost(  ) :
    m_a(NULL)
    , m_dest(NULL)
    , m_size(0)
    , m_capacity(0)
    , m_maxBlocks(0)
    , m_maxThreadsPerBlock(0)
    , m_status(true)
    , m_gen()
{
    init();
}

void SquareHost::init( ) {
    curandStatus_t err =  curandCreateGenerator( &m_gen, CURAND_RNG_PSEUDO_MTGP32 );
    if( err != CURAND_STATUS_SUCCESS ) {
        std::cerr << "ERROR: Failed to create generator: " << err << std::endl;
        m_status = false;
        return;
    }
}

void SquareHost::operator()( unsigned int s, seed_type seed ) {
    resize( s );

    if( !good() ) return;

    if( curandSetPseudoRandomGeneratorSeed( m_gen, seed ) != CURAND_STATUS_SUCCESS ) {
        m_status = false;
        return;
    }

    curandGenerate( m_gen, m_dest, s );

    unsigned int bcount = (s / 1024);
    if( s % 1024 ) { ++bcount; }

    dim3 grid(1, bcount, 1), block(256, 4, 1);

    computeHW<<< grid, block >>>( m_dest, s );

    cudaMemcpy( m_a, m_dest, s * sizeof(int_type), cudaMemcpyDeviceToHost); // blocks host until copy is complete
}

void SquareHost::resize( unsigned int s ) {
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
            m_status = false;
        }

        m_capacity = s;
    }

    m_size = s;
}

bool SquareHost::good() const { return m_status; }


SquareHost::~SquareHost() {
    if( m_status ) { curandDestroyGenerator( m_gen ); }

    if( m_a ) free(m_a);

    if( m_dest ) cudaFree(m_dest);
}

std::ostream & operator<<( std::ostream & out, const SquareHost & rhs ) {
    if( !rhs.good() ) {
        out << "BAD STATE";
    } else {
        if( rhs.m_size ) {
            unsigned int i = 0;
            out << rhs.m_a[i++];
            while( i < rhs.m_size ) {
                out << "," << rhs.m_a[i++];
            }
        }
    }

    return out;
}
