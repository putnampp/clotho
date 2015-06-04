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
#include "cuda_mt19937.h"

#include <curand_mtgp32_host.h>
#include <curand_mtgp32dc_p_11213.h>

#include <iostream>

dim3 cuda_mt19937::makeBlockDimension( unsigned int max ) {
    dim3 res( THREADS_PER_STATE, max / THREADS_PER_STATE, 1 );

    if( max % THREADS_PER_STATE ) {
        ++res.y;
    }
    return res;
}

cuda_mt19937::cuda_mt19937( boost::random::mt19937 * gen, unsigned int s ) :
    m_states( NULL )
    , m_size( 0 )
    , m_capacity(0)
    , m_status( true )
    , m_dKernelParams( NULL )
    , m_rng( gen ) 
{
    resize(s);
}

cuda_mt19937::pointer cuda_mt19937::getStates() {
    return m_states;
}

bool cuda_mt19937::good() const {
    return m_status;
}

void cuda_mt19937::updateStates( unsigned int s ) {
    if( !resize( s ) ) return;

    pointer state = m_states;
    while( s ) {
        unsigned int nStates = (( s >= 200 ) ? 200 : s );

        if( curandMakeMTGP32KernelState( state, mtgp32dc_params_fast_11213, m_dKernelParams, nStates, m_dist( *m_rng ) ) != CURAND_STATUS_SUCCESS ) {
            std::cerr << "ERROR: Unable to make state (" << m_states << ", " << state << ", " << sizeof(state_t) << ")" << std::endl;
            m_status = false;
            return;
        }

        state += nStates;
        s -= nStates;
    }
}

bool cuda_mt19937::resize( unsigned int s ) {
    reserve( s );
    m_size = s;
    
    return m_status;
}

bool cuda_mt19937::reserve( unsigned int s ) {
    if( !good() ) return m_status;

    if( s > m_capacity ) {
//        std::cerr << "Increasing state capacity " << m_capacity << " -> " << s << std::endl;
        if( m_capacity > 0 ) {
            size_t tmp_size = s * sizeof( state_t );
            pointer t = m_states;
            if( cudaMalloc( (void **) &m_states, tmp_size) != cudaSuccess ) {
                std::cerr << "ERROR: Unable to allocate state memory: " << tmp_size << std::endl;
                m_status = false;
                return m_status;
            }
            tmp_size = m_capacity * sizeof( state_t );
            if( cudaMemcpy( m_states, t, tmp_size, cudaMemcpyDeviceToDevice ) != cudaSuccess ) {
                std::cerr << "ERROR: Unable to copy previous states: " << tmp_size << std::endl;
                m_status = false;
                return m_status;
            }
            cudaFree( t );

//            tmp_size = s * sizeof(mtgp32_kernel_params);
//            mtgp32_kernel_params * c = m_dKernelParams;
//            if( cudaMalloc( (void **) &m_dKernelParams, tmp_size ) != cudaSuccess ) {
//                std::cerr << "ERROR: Unable allocate new parameters (" << tmp_size << " bytes)" << std::endl;
//                m_status = false;
//                return;
//            }
//
//            tmp_size = m_capacity * sizeof( mtgp32_kernel_params );
//            if( cudaMemcpy( m_dKernelParams, c, tmp_size, cudaMemcpyDeviceToDevice ) != cudaSuccess ) {
//                std::cerr << "ERROR: Unable copy parameters (" << tmp_size << " bytes)" << std::endl;
//                m_status = false;
//                return; 
//            }
//            cudaFree( c );
        } else {
            size_t tmp_size = s * sizeof( state_t );
            if( cudaMalloc( (void **) &m_states, tmp_size ) != cudaSuccess ) {
                std::cerr << "ERROR: Unable to allocate states (" << tmp_size << " bytes)" << std::endl;
                m_status = false;
                return m_status;
            }

//            tmp_size = s * sizeof(mtgp32_kernel_params);
//
//            if( cudaMalloc( (void **) &m_dKernelParams, tmp_size ) != cudaSuccess ) {
//                std::cerr << "ERROR: Unable to allocate kernel parameters (" << tmp_size << " bytes) " << std::endl;
//                m_status = false;
//                return;
//            }
        }

        if( m_dKernelParams == NULL ) {
            if( cudaMalloc( (void **)&m_dKernelParams, sizeof( mtgp32_kernel_params) ) != cudaSuccess ) {
                std::cerr << "Unable to allocate kernel_params (" << sizeof( mtgp32_kernel_params ) << " bytes)" << std::endl;
                m_status = false;
                return m_status;
            }

            if( curandMakeMTGP32Constants( mtgp32dc_params_fast_11213, m_dKernelParams ) != CURAND_STATUS_SUCCESS ) {
                std::cerr << "ERROR: Unable to make constants (" << m_dKernelParams << ")" << std::endl;
                m_status = false;
                return m_status;
            }
        }
        
        /*pointer state = m_states + m_capacity;
        while( m_capacity < s ) {
            unsigned int nStates = s - m_capacity;
            if( nStates > 200 ) {
                nStates = 200;
            }
            if( curandMakeMTGP32KernelState( state, mtgp32dc_params_fast_11213, m_dKernelParams, nStates, m_dist( *m_rng ) ) != CURAND_STATUS_SUCCESS ) {
                std::cerr << "ERROR: Unable to make state (" << m_states << ", " << state << ", " << sizeof(state_t) << ")" << std::endl;
                m_status = false;
                return;
            }

            state += nStates;
            m_capacity += nStates;
        }*/
        m_capacity = s;
    }
    return m_status;
}

cuda_mt19937::~cuda_mt19937() {
    if( m_states ) cudaFree( m_states );
    if( m_dKernelParams ) cudaFree( m_dKernelParams );
}
