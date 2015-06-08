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
#ifndef SQUARE_CUDA_STREAM_H_
#define SQUARE_CUDA_STREAM_H_

#include <ostream>

#include <cuda.h>

#include <boost/random/mersenne_twister.hpp>

class SquareStream {
public:
    typedef unsigned long int_type;

    static const unsigned int BYTES_PER_STREAM = (1 << 22); // 22 => 4MB/Stream
    static const unsigned int INT_PER_STREAM = BYTES_PER_STREAM / (2 * sizeof( int_type ));
    static const unsigned int MAX_THREADS = 1024;

    SquareStream( boost::random::mt19937 & rng, unsigned int nStreams = 16);

    size_t size() const;

    void operator()( unsigned int s );

//    void random_list();

    bool good() const;

    friend std::ostream & operator<<( std::ostream &, const SquareStream & rhs );

    virtual ~SquareStream();
protected:

    void init();

    void method1( unsigned int s );
    void method2( unsigned int s );

/**
 * Resize both host and device memory arrays to
 * specific size.
 *
 * If size is greater than current capacity
 * then both arrays are reallocated with size as new capacity
 * and all existing data is lost
 */
    void resize( unsigned int s );

    int_type *  m_hostMem, * m_devMem;
    cudaStream_t * m_streams;
    int_type ** m_hBuffer, ** m_dBuffer;
    int_type * m_streamSizes;
    size_t      m_size, m_capacity;
    unsigned int m_nStreams;
    int m_maxBlocks, m_maxThreadsPerBlock;
    bool    m_status;

    boost::random::mt19937 * m_rng;
};

#endif  // SQUARE_CUDA_STREAM_H_
