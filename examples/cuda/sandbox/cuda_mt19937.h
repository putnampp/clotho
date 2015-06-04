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
#ifndef CUDA_MT19937_H_
#define CUDA_MT19937_H_

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include <cuda.h>
#include <curand_kernel.h>

class cuda_mt19937 {
public:
    typedef curandStateMtgp32_t  state_t;
    typedef state_t * pointer;

    static const unsigned int THREADS_PER_STATE = 256;

    static dim3 makeBlockDimension( unsigned int max_threads );

    cuda_mt19937( boost::random::mt19937 * gen, unsigned int s = 0 );

    pointer getStates();

    bool good() const;

    void updateStates( unsigned int s );

    virtual ~cuda_mt19937();
protected:
    bool resize( unsigned int s );
    bool reserve( unsigned int s );

    pointer m_states;
    size_t  m_size, m_capacity;
    bool m_status;

    mtgp32_kernel_params * m_dKernelParams;

    boost::random::mt19937 * m_rng;
    boost::random::uniform_int_distribution< unsigned long > m_dist;
};
#endif  // CUDA_MT19937_H_
