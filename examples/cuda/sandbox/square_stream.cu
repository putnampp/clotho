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
#include "square_stream.h"

#include <stdio.h>
#include <iostream>
#include <cassert>

#include "popcount_kernel.h"

__global__ void square( SquareStream::int_type * a, unsigned int N ) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    //__shared__ SquareStream::int_type buffer[ SquareStream::MAX_THREADS ];
    
    //buffer[ threadIdx.x ] = a[idx]; // assumes the device memory is padded
    //buffer[threadIdx.x] *= buffer[threadIdx.x];
    //a[idx] = buffer[threadIdx.x];
    a[idx] *= a[idx];
}
/*
template < class I >
__device__ I popcountGPU( I a );

template <>
__device__ unsigned long popcountGPU( unsigned long x ) {
    x -= (x >> 1) & m1;
    x = (x & m2) + ((x >> 2) & m2);
    x = (x + (x >> 4) ) & m4;
    return ( x * h01) >> 56;
}

template <>
__device__ unsigned int popcountGPU( unsigned int x ) {
    x -= (x >> 1) & (unsigned int) m1;
    x = (x & (unsigned int) m2) + ((x >> 2) & (unsigned int) m2);
    x = (x + (x >> 4) ) & (unsigned int) m4;
    return ( x * (unsigned int) h01) >> 24;
}

__global__ void computeHW( SquareStream::int_type * a, unsigned int N ) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    a[idx] = popcountGPU(a[idx]);
}*/

SquareStream::SquareStream( boost::random::mt19937 & rng, unsigned int nStreams ) :
    m_hostMem(NULL)
    , m_devMem(NULL)
    , m_streams(NULL)
    , m_hBuffer(NULL)
    , m_dBuffer(NULL)
    , m_streamSizes(NULL)
    , m_size(0)
    , m_capacity(0)
    , m_nStreams( nStreams )
    , m_maxBlocks(0)
    , m_maxThreadsPerBlock(0)
    , m_status(true)
    , m_rng( &rng )
{
    init();
}

SquareStream::SquareStream( boost::random::mt19937 * rng, unsigned int nStreams ) :
    m_hostMem(NULL)
    , m_devMem(NULL)
    , m_streams(NULL)
    , m_hBuffer(NULL)
    , m_dBuffer(NULL)
    , m_streamSizes(NULL)
    , m_size(0)
    , m_capacity(0)
    , m_nStreams( nStreams )
    , m_maxBlocks(0)
    , m_maxThreadsPerBlock(0)
    , m_status(true)
    , m_rng( rng )
{
    init();
}
void SquareStream::init() {
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

    if( m_nStreams ) {
        m_streams = (cudaStream_t *) malloc( m_nStreams * sizeof(cudaStream_t) );
        cudaStream_t * tmp = m_streams;
        unsigned int i = m_nStreams;
        while( i-- ) {
            cudaStreamCreateWithFlags( tmp, cudaStreamNonBlocking );
            ++tmp;
        }
        
        i = m_nStreams * BYTES_PER_STREAM;

        if( cudaMalloc( (void **) &m_devMem, i ) != cudaSuccess ) {
            std::cerr << "Unable to allocate device memory: " << std::endl;
            m_status = false;
            return;
        }

        m_hBuffer = (int_type **)malloc( m_nStreams * sizeof(int_type*));
        m_dBuffer = (int_type **)malloc( m_nStreams * sizeof(int_type*));
        m_streamSizes = (int_type *)malloc( m_nStreams * sizeof( int_type) );
    }
}

SquareStream::~SquareStream() {
    if( m_hostMem ) free(m_hostMem);
    if( m_devMem ) cudaFree( m_devMem );
    if( m_streams ) {
        while( m_nStreams-- ) {
            cudaStreamDestroy( m_streams[m_nStreams] );
        }
        free(m_streams);
    }
    if( m_hBuffer ) free( m_hBuffer );
    if( m_dBuffer ) free( m_dBuffer );
    if( m_streamSizes ) free( m_streamSizes );
}

size_t SquareStream::size() const { return m_size; }

bool SquareStream::good() const {
    return m_status;
}

void SquareStream::operator()( unsigned int s ) {
    resize( s );

    assert(good());

    method2( s );
}

void SquareStream::method1( unsigned int s ) {
    unsigned int bcount = INT_PER_STREAM / MAX_THREADS; // Block/Stream
    if( INT_PER_STREAM % MAX_THREADS ) { ++bcount; }
    
    int_type * hMem = m_hostMem;
    do {
        int_type * dMem = m_devMem;

        unsigned int used_streams = 0;
        while( s && used_streams < m_nStreams ) {
            m_hBuffer[ used_streams ] = hMem;
            m_dBuffer[ used_streams ] = dMem;

            unsigned int n = ((s >= INT_PER_STREAM ) ? INT_PER_STREAM : s);
            m_streamSizes[ used_streams ] = n;

            s -= n;
            dMem += n;
            while( n-- ) {
                (*hMem++) = (*m_rng)();
            }

            cudaMemcpyAsync( m_dBuffer[used_streams], m_hBuffer[used_streams], m_streamSizes[used_streams] * sizeof( int_type ), cudaMemcpyHostToDevice, m_streams[used_streams] );

            ++used_streams;
        }

        unsigned int i = 0;
        while( i < used_streams ) {
            computeHW<<< bcount, MAX_THREADS, 0, m_streams[i]>>>( m_dBuffer[i], m_streamSizes[i] );
            ++i;
        }

        while( used_streams-- ) {
            cudaMemcpyAsync( m_hBuffer[used_streams], m_dBuffer[used_streams], m_streamSizes[ used_streams ] * sizeof(int_type), cudaMemcpyDeviceToHost, m_streams[used_streams]);
       }
    } while ( s );
    cudaDeviceSynchronize();
}

void SquareStream::method2( unsigned int s ) {
    unsigned int bcount = INT_PER_STREAM / MAX_THREADS; // Block/Stream
    if( INT_PER_STREAM % MAX_THREADS ) { ++bcount; }
    
    int_type * hMem = m_hostMem;
    do {
        int_type * dMem = m_devMem;

        unsigned int used_streams = 0;
        while( s && used_streams < m_nStreams ) {
            m_hBuffer[ used_streams ] = hMem;
            m_dBuffer[ used_streams ] = dMem;

            unsigned int n = ((s >= INT_PER_STREAM ) ? INT_PER_STREAM : s);
            m_streamSizes[ used_streams ] = n;

            s -= n;
            dMem += n;
            while( n-- ) {
                (*hMem++) = (*m_rng)();
            }
            ++used_streams;
        }

        unsigned int i = 0;
        while( i < used_streams ) {
            cudaMemcpyAsync( m_dBuffer[i], m_hBuffer[i], m_streamSizes[i] * sizeof(int_type), cudaMemcpyHostToDevice, m_streams[i] );
            computeHW<<< bcount, MAX_THREADS,0, m_streams[i]>>>( m_dBuffer[i], m_streamSizes[i] );

            cudaMemcpyAsync( m_hBuffer[i], m_dBuffer[i], m_streamSizes[i] * sizeof(int_type), cudaMemcpyDeviceToHost, m_streams[i]);
            ++i;
        }

    } while ( s );
}

void SquareStream::resize( unsigned int s ) {
    if( !good() ) return;

    if ( s > m_capacity ) {
        if( m_hostMem ) {
            free( m_hostMem );
        }

        size_t byte_size = s * sizeof(int_type );
        m_hostMem = (int_type *) malloc( byte_size );

        m_capacity = s;
    }

    m_size = s;
}

std::ostream & operator<<( std::ostream & out, const SquareStream & rhs ) {
    if( rhs.good() ) {
        for( unsigned int i = 0; i < rhs.size(); ++i ) {
            out << i << " -> " << rhs.m_hostMem[i] << "\n";
        }
    } else {
        out << "BAD STATE" << "\n";
    }
    return out;
}
