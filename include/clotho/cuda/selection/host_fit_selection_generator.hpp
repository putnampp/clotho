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
#ifndef HOST_FIT_SELECTION_GENERATOR_HPP_
#define HOST_FIT_SELECTION_GENERATOR_HPP_

#include "clotho/recombination/sequence_bias_parameter.hpp"

#include <boost/random/discrete_distribution.hpp>
#include <boost/random/bernoulli_distribution.hpp>

class HostSelectionGenerator : public clotho::utility::iStateObject {
public:

    typedef sequence_bias_parameter< double > sequence_bias_type;

    HostSelectionGenerator( boost::property_tree::ptree & config ) :
        m_seq_bias( config )
        , m_hSelection(NULL)
        , m_dSelection(NULL)
        , m_size(0)
        , m_capacity(0)
        , m_dSize(0)
        , m_dCapacity(0)
    {
        cudaStreamCreate( &m_upload );
    }

    void get_state( boost::property_tree::ptree & s ) {
        boost::property_tree::ptree d;

        for( unsigned int i = 0; i < m_size; ++i ) {
            clotho::utility::add_value_array( d, m_hSelection[ i ] );
        }
        s.put( "size", m_size );
        s.put( "capacity", m_capacity );
        s.put_child( "data", d );
    }

    template < class RNG, class RealType, class IntType >
    void operator()( RNG & rng, HostPopulationSpace< RealType, IntType > * parent, HostPopulationSpace< RealType, IntType > * offspring ) {

        typedef HostPopulationSpace< RealType, IntType > population_space_type;
        typedef typename population_space_type::fitness_type fitness_type;

        resize( offspring->getSequenceCount() );

        if( parent->getIndividualCount() == 0 ) {
            memset( m_hSelection, 0, m_size * sizeof(unsigned int));
        } else {
            boost::random::discrete_distribution< unsigned int, fitness_type > disc_dist( parent->getHostFitness(), parent->getHostFitness() + parent->getIndividualCount() );
            boost::random::bernoulli_distribution< double > bern( m_seq_bias.m_bias );

            for( unsigned int i = 0; i < m_size; ++i ) {
                m_hSelection[ i ] = 2 * disc_dist( rng );
                if( bern( rng ) ) {
                    m_hSelection[ i ] += 1;
                }
            }
        }
    }

    unsigned int * getDeviceList() {
        return m_dSelection;
    }

    size_t getDeviceSize() const {
        return m_dSize;
    }

    void resizeDevice() {
        if( m_size > m_dCapacity ) {
            if( m_dSelection != NULL ) {
                assert( cudaFree( m_dSelection ) == cudaSuccess);
                m_dSelection = NULL;
            }

            assert( cudaMalloc( (void **) & m_dSelection, m_size * sizeof( unsigned int ) ) == cudaSuccess );
            m_dCapacity = m_size;
        }
        m_dSize = m_size;
    }

    void updateDevice( ) {
        resizeDevice();
        if( m_dSize > 0 ) {
            assert( m_dSelection != NULL );
            assert( m_hSelection != NULL );

            cudaError_t err = cudaMemcpy( m_dSelection, m_hSelection, m_dSize * sizeof(unsigned int), cudaMemcpyHostToDevice);
            if( err != cudaSuccess ) {
                std::cerr << "ERROR: " << cudaGetErrorString( err ) << std::endl;
                assert( false );
            }
        }
    }

    void updateDeviceAsync( ) {
        resizeDevice();
        if( m_dSize > 0 ) {
            assert( m_dSelection != NULL );
            assert( m_hSelection != NULL );

            cudaError_t err = cudaMemcpyAsync( m_dSelection, m_hSelection, m_dSize * sizeof(unsigned int), cudaMemcpyHostToDevice, m_upload);
            if( err != cudaSuccess ) {
                std::cerr << "ERROR: " << cudaGetErrorString( err ) << std::endl;
                assert( false );
            }
        }
    }

    void sync() {
        cudaStreamSynchronize( m_upload );
    }

    virtual ~HostSelectionGenerator() {
        if( m_hSelection != NULL ) {
            delete [] m_hSelection;
        }
        if( m_dSelection != NULL ) {
            cudaFree( m_dSelection );
        }

        cudaStreamDestroy( m_upload );
    }

protected:

    void resize( size_t N ) {
        if( N > m_capacity ) {
//            std::cerr << "Fitness Resizing: " << N << std::endl;
            if( m_hSelection != NULL ) {
                delete [] m_hSelection;
            }

            m_hSelection = new unsigned int[ N ];
            m_capacity = N;
        }
        m_size = N;
    }

    sequence_bias_type m_seq_bias;
    unsigned int * m_hSelection, * m_dSelection;

    size_t m_size, m_capacity;
    size_t m_dSize, m_dCapacity;

    cudaStream_t m_upload;
};

#endif  // HOST_FIT_SELECTION_GENERATOR_HPP_
