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
#ifndef HOST_SEQUENCE__WEIGHT_HPP_
#define HOST_SEQUENCE__WEIGHT_HPP_

#include "clotho/cuda/analysis/sequence_weight_kernel.hpp"

class HostSequenceWeight : public clotho::utility::iStateObject {
public:

    HostSequenceWeight() :
        m_hSeqWeights( NULL )
        , m_dSeqWeights( NULL )
        , m_size(0)
        , m_capacity(0)
    {}

    void get_state( boost::property_tree::ptree & s ) {
        updateHost();

        boost::property_tree::ptree d;
        for( unsigned int i = 0; i < m_size; ++i ) {
            clotho::utility::add_value_array( d, m_hSeqWeights[ i ] );
        }

        s.put("size", m_size);
        s.put("capacity", m_capacity );
        s.add_child("data", d);
    }

    void updateHost() {
        assert( cudaMemcpy( m_hSeqWeights, m_dSeqWeights, m_size * sizeof( unsigned int ), cudaMemcpyDeviceToHost ) == cudaSuccess );
    }

    template < class RealType, class IntType >
    void operator()( HostPopulationSpace< RealType, IntType > * pop ) {
        resize( pop->getSequenceCount() );

        dim3 blocks( pop->getSequenceCount(), 1, 1), threads( 32, 32, 1 );

        evaluate_sequence_weights<<< blocks, threads >>>( pop->getDeviceSequences(), m_dSeqWeights, pop->getBlocksPerSequence());
    }

    virtual ~HostSequenceWeight() {
        if( m_hSeqWeights != NULL ) {
            delete [] m_hSeqWeights;
            cudaFree( m_dSeqWeights );
        }
    }

protected:
    void resize( unsigned int N ) {
        if( N > m_capacity ) {
            if( m_hSeqWeights != NULL ) {
                delete [] m_hSeqWeights;
                cudaFree( m_dSeqWeights );
            }
            assert( cudaMalloc( (void **) &m_dSeqWeights, N * sizeof( unsigned int ) ) == cudaSuccess );

            m_hSeqWeights = new unsigned int[ N ];

            m_capacity = N;
        }
        m_size = N;
    }

    unsigned int * m_hSeqWeights, * m_dSeqWeights;
    unsigned int m_size, m_capacity;
};

#endif  // HOST_SEQUENCE__WEIGHT_HPP_
