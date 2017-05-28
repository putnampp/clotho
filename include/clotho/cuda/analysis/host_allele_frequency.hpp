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
#ifndef HOST_ALLELE_FREQUENCY_HPP_
#define HOST_ALLELE_FREQUENCY_HPP_

#include "clotho/cuda/analysis/allele_frequency_kernels.hpp"

class HostAlleleFrequency : public clotho::utility::iStateObject {
public:

    HostAlleleFrequency( ):
        m_hAllFreq( NULL )
        , m_dAllFreq( NULL )
        , m_size(0)
        , m_capacity(0)
    {}

    void get_state( boost::property_tree::ptree & s ) {
        updateHost();
    }

    template < class RealType, class IntType >
    void operator()( HostPopulationSpace< RealType, IntType > * pop, unsigned int N ) {
        resize( N );

        unsigned int b_x = N / 1024 + ((N % 1024) ? 1 : 0);

        dim3 blocks( b_x, 1, 1), threads( 32, 32, 1 );

        count_alleles<<< blocks, threads >>>( pop->getDeviceSequences(), m_dAllFreq, pop->getSequenceCount(), pop->getBlocksPerSequence(), N );
    }

    void updateHost() {
        assert( cudaMemcpy( m_hAllFreq, m_dAllFreq, m_size * sizeof(unsigned int), cudaMemcpyDeviceToHost) == cudaSuccess );
    }

    virtual ~HostAlleleFrequency() {
        if( m_hAllFreq != NULL ) {
            delete [] m_hAllFreq;
            cudaFree( m_dAllFreq );
        }
    }

protected:
    void resize( unsigned int N ) {
        if( N > m_capacity ) {
            if( m_hAllFreq != NULL ) {
                delete [] m_hAllFreq;
                cudaFree( m_dAllFreq );
            }

            assert( cudaMalloc( (void **) &m_dAllFreq, N * sizeof( unsigned int) ) == cudaSuccess );
            m_hAllFreq = new unsigned int [ N * sizeof( unsigned int ) ];

            m_capacity = N;
        }
        m_size = N;
    }

    unsigned int * m_hAllFreq, * m_dAllFreq;

    unsigned int m_size, m_capacity;
};

#endif  // HOST_ALLELE_FREQUENCY_HPP_
