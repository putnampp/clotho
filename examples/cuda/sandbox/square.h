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
#ifndef SQUARE_CUDA_H_
#define SQUARE_CUDA_H_

#include <ostream>

#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>

#include <curand_mtgp32.h>

#include <boost/random/mersenne_twister.hpp>
#include "cuda_mt19937.h"

class Square {
public:
    typedef unsigned long int_type;
    typedef unsigned long seed_type;

    typedef cuda_mt19937 curand_state_type;

    Square( boost::random::mt19937 & rng);

    size_t size() const;

    void operator()( unsigned int s );

//    void random_list();

    bool good() const;

    friend std::ostream & operator<<( std::ostream &, const Square & rhs );

    virtual ~Square();
protected:

    void init();

/**
 * Resize both host and device memory arrays to
 * specific size.
 *
 * If size is greater than current capacity
 * then both arrays are reallocated with size as new capacity
 * and all existing data is lost
 */
    void resize( unsigned int s );

    void initializePRNG( curandStateMtgp32_t * );

    int_type *  m_a, * m_dest;
    size_t      m_size, m_capacity;
    int m_maxBlocks, m_maxThreadsPerBlock;
    bool    m_status;

    curand_state_type m_dStates;
};

#endif  // SQUARE_CUDA_H_
