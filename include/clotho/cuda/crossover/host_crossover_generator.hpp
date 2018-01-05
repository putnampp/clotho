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
#ifndef HOST_CROSSOVER_GENERATOR_HPP_
#define HOST_CROSSOVER_GENERATOR_HPP_

#include "clotho/recombination/recombination_rate_parameter.hpp"

#include <boost/random/poisson_distribution.hpp>
#include <boost/random/uniform_01.hpp>

#include "clotho/cuda/crossover/crossover_kernel_build_mask_impl.hpp"
#include "clotho/cuda/crossover/crossover_kernel.hpp"

#include "clotho/cuda/data_spaces/sequence_space/sequence_kernels.hpp"

#include "clotho/utility/algorithm_version.hpp"

template < class EventType >
class HostCrossoverGenerator {
public:

    static const unsigned int MAX_MUTATIONS_PER_SEQUENCE = 32;

    typedef EventType   event_type;
    
    typedef boost::random::poisson_distribution< unsigned int, double > event_distribution_type;
    typedef boost::random::uniform_01< event_type > location_distribution_type;

    typedef typename build_crossover_mask::block_sequence_kernel build_mask_version;

    HostCrossoverGenerator( boost::property_tree::ptree & config ) :
        m_recomb( config )
        , m_hEventPool(NULL)
        , m_dEventPool(NULL)
        , m_hEventDist(NULL)
        , m_dEventDist(NULL)
        , m_pool_offset(0)
        , m_pool_capacity(0)
        , m_dist_size(0)
        , m_dist_capacity(0)
        , m_event_pool_reset( false )
    {
        cudaStreamCreate( &m_maskStream );
    }

    template < class RNG, class RealType, class IntType >
    void initialize( RNG & rng, HostPopulationSpace< RealType, IntType > * pop ) {
        generateEvents( rng, pop->getSequenceCount() );
    }

    template < class RealType, class IntType >
    void buildRecombinationMask( HostAlleleSpace< RealType > & alleles, HostPopulationSpace< RealType, IntType > * offspring ) {
       updateDevice();

        if( alleles.getDeviceMaxAlleles() == 0 ) {
//            std::cerr << "No device alleles" << std::endl;

            clear_sequence_space::execute( offspring->getDeviceSequences(), offspring->getSequenceCount(), offspring->getBlocksPerSequence() );

        } else {
            build_mask_version v;
            build_crossover_mask::execute( alleles.getDeviceLocations(), m_dEventPool, m_dEventDist, offspring->getDeviceSequences(), offspring->getSequenceCount(), offspring->getBlocksPerSequence(), alleles.getDeviceAlleleCount(), &v );
        }
    }

    template < class RealType, class IntType >
    void buildRecombinationMaskAsync( HostAlleleSpace< RealType > & alleles, HostPopulationSpace< RealType, IntType > * offspring ) {
       updateDeviceAsync();

        if( alleles.getDeviceMaxAlleles() == 0 ) {
            clear_sequence_space::execute( offspring->getDeviceSequences(), offspring->getSequenceCount(), offspring->getBlocksPerSequence(), m_maskStream );
        } else {
            build_mask_version v;
            build_crossover_mask::execute( alleles.getDeviceLocations(), m_dEventPool, m_dEventDist, offspring->getDeviceSequences(), offspring->getSequenceCount(), offspring->getBlocksPerSequence(), alleles.getDeviceAlleleCount(), m_maskStream, &v );
        }
    }

    template < class RealType, class IntType >
    void performCrossover( HostPopulationSpace< RealType, IntType > * parents, HostPopulationSpace< RealType, IntType > * offspring, HostSelectionGenerator & sel ) {
        if( parents->getDeviceSequences() == NULL ) {
            // some what redundant as the only time the parents will be null
            // is when the there are no allele locations
            // so the buildMask function should have already cleared this space
            // this just adds to the start-up cost
            clear_sequence_space::execute( offspring->getDeviceSequences(), offspring->getSequenceCount(), offspring->getBlocksPerSequence() );
        } else {
            assert( sel.getDeviceSize() == offspring->getSequenceCount() );

            crossover::execute( parents->getDeviceSequences(), offspring->getDeviceSequences(), sel.getDeviceList(), offspring->getSequenceCount(), parents->getBlocksPerSequence(), offspring->getBlocksPerSequence() );
        }
    }

    void updateDevice() {
        if( m_event_pool_reset ) {
            assert( cudaMemcpy( m_dEventPool, m_hEventPool, m_pool_capacity * sizeof( event_type ), cudaMemcpyHostToDevice ) == cudaSuccess );
            m_event_pool_reset = false;
        }
        assert( cudaMemcpy( m_dEventDist, m_hEventDist, m_dist_size * sizeof( unsigned int ), cudaMemcpyHostToDevice ) == cudaSuccess );
    }

    void updateDeviceAsync() {
        if( m_event_pool_reset ) {
            assert( cudaMemcpyAsync( m_dEventPool, m_hEventPool, m_pool_capacity * sizeof( event_type ), cudaMemcpyHostToDevice, m_maskStream ) == cudaSuccess );
            m_event_pool_reset = false;
        }
        assert( cudaMemcpyAsync( m_dEventDist, m_hEventDist, m_dist_size * sizeof( unsigned int ), cudaMemcpyHostToDevice, m_maskStream ) == cudaSuccess );
    }

    void sync() {
        assert( cudaStreamSynchronize( m_maskStream ) == cudaSuccess );
    }

    virtual ~HostCrossoverGenerator() {
        if( m_hEventDist != NULL ) {
            delete [] m_hEventDist;
            cudaFree( m_dEventDist );

            delete [] m_hEventPool;
            cudaFree( m_dEventPool );
        }

        cudaStreamDestroy( m_maskStream );
    }

protected:

    template < class RNG >
    void generateEventPool( RNG & rng ) {
        location_distribution_type loc_dist;

        for( unsigned int i = 0; i < m_pool_capacity; ++i ) {
            m_hEventPool[ i ] = loc_dist( rng );
        }
    }

    template < class RNG >
    void generateEvents( RNG & rng, unsigned int N ) {
        resizeDistribution( N );
        assert( N < m_dist_size );

        event_distribution_type events( m_recomb.m_rho );

        m_hEventDist[ 0 ] = m_pool_offset;
        for( unsigned int i = 1; i <= N; ++i ) {
            unsigned int E = events( rng );

            m_pool_offset += E;

            m_hEventDist[ i ] = m_pool_offset;  
        }

        if( m_pool_offset >= m_pool_capacity ) {
            // just in case sequence space has increased since previous call
            resizePool( N );

            assert( m_dEventPool != NULL );
            assert( (m_pool_offset - m_hEventDist[ 0 ]) < m_pool_capacity );

            generateEventPool( rng );

            unsigned int sub_val = m_hEventDist[ 0 ];
            for( unsigned int i = 0; i <= N; i++ ) {
                m_hEventDist[ i ] -= sub_val;
            }

            m_pool_offset = m_hEventDist[ N ];
            m_event_pool_reset = true;
        }
    }

    void resizeDistribution( unsigned int N ) {
        if( N > m_dist_capacity ) {
            if( m_dEventDist != NULL ) {
                delete [] m_hEventDist;
                cudaFree( m_dEventDist );
            }

            m_hEventDist = new unsigned int[ N + 1 ];
            assert( cudaMalloc( (void **) &m_dEventDist, (N + 1) * sizeof( unsigned int) ) == cudaSuccess );
            m_dist_capacity = N + 1;
        }
        m_dist_size = N + 1;
    }

    void resizePool( unsigned int N ) {
        // always resize relative to the population size
        unsigned int pool_cap = MAX_MUTATIONS_PER_SEQUENCE * N;
        if( pool_cap > m_pool_capacity ) {
            if( m_hEventPool != NULL ) {
                delete [] m_hEventPool;
                cudaFree( m_dEventPool );
            }

            m_hEventPool = new event_type[ pool_cap ];
            assert( cudaMalloc( (void **) &m_dEventPool, pool_cap * sizeof( event_type ) ) == cudaSuccess );

            m_pool_capacity = pool_cap;
        }
    }

    recombination_rate_parameter< double > m_recomb;

    event_type * m_hEventPool, * m_dEventPool;
    unsigned int * m_hEventDist, * m_dEventDist;

    unsigned int m_pool_offset, m_pool_capacity;
    unsigned int m_dist_size, m_dist_capacity;

    cudaStream_t m_maskStream;
    bool m_event_pool_reset;
};

#endif  // HOST_CROSSOVER_GENERATOR_HPP_
