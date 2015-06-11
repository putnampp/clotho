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
#ifndef SQUARE_HOST_CUDA_H_
#define SQUARE_HOST_CUDA_H_

#include <ostream>

#include <cuda.h>
#include <curand.h>

class SquareHost {
public:
    typedef unsigned int int_type;
    typedef unsigned long long seed_type;

    SquareHost(  );

    size_t size() const;

    void operator()( unsigned int s, seed_type seed );

    bool good() const;

    friend std::ostream & operator<<( std::ostream &, const SquareHost & rhs );

    virtual ~SquareHost();
protected:

    void init( );

/**
 * Resize both host and device memory arrays to
 * specific size.
 *
 * If size is greater than current capacity
 * then both arrays are reallocated with size as new capacity
 * and all existing data is lost
 */
    void resize( unsigned int s );

    int_type *  m_a, * m_dest;
    size_t      m_size, m_capacity;
    int m_maxBlocks, m_maxThreadsPerBlock;
    bool    m_status;

    curandGenerator_t m_gen;
};

#endif  // SQUARE_HOST_CUDA_H_
