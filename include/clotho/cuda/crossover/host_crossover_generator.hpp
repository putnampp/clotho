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

template < class EventType >
class HostCrossoverGenerator {
public:

    typedef EventType   event_type;
    
    typedef boost::random::poisson_distribution< unsigned int, double > event_distribution_type;
    typedef boost::random::uniform_01< event_type > location_distribution_type;

    HostCrossoverGenerator( boost::property_tree::ptree & config ) :
        m_recomb( config )
        , m_hEventPool(NULL)
        , m_dEventPool(NULL)
        , m_hEventDist(NULL)
        , m_dEventDist(NULL)
        , m_pool_size(0)
        , m_pool_capacity(0)
        , m_dist_size(0)
        , m_dist_capacity(0)
    {}

    template < class RNG, class RealType, class IntType >
    void initialize( RNG & rng, HostPopulationSpace< RealType, IntType > * pop ) {
        location_distribution_type loc_dist;
        event_distribution_type events( m_recomb.m_rho );

        resize( pop->getSequenceCount() );

        m_hEventDist[ 0 ] = 0;
        for( unsigned int i = 0; i < m_dist_size; ++i ) {
            unsigned int N = events( rng );
            while( N > 0 ) {
                m_hEventPool[ m_pool_size++ ] = loc_dist( rng );
                --N;
            }

            m_hEventDist[ i + 1 ] = m_pool_size;  
        }
        updateDevice();
    }

    template < class RealType, class IntType >
    void buildRecombinationMask( HostAlleleSpace< RealType > & alleles, HostPopulationSpace< RealType, IntType > * offspring ) {

        dim3 blocks( 1, 1, 1), threads( 1,1,1);

        computeBlockThreadDims( offspring->getSequenceCount(), offspring->getBlocksPerSequence(), blocks, threads );

        build_crossover_mask<<< blocks, threads >>>( alleles.getDeviceLocations(), m_dEventPool, m_dEventDist, offspring->getDeviceSequences(), offspring->getBlocksPerSequence(), alleles.getAlleleCount() );
        CHECK_LAST_KERNEL_EXEC
    }

    template < class RealType, class IntType >
    void performCrossover( HostPopulationSpace< RealType, IntType > * parents, HostPopulationSpace< RealType, IntType > * offspring, HostSelectionGenerator & sel ) {

        sel.updateDevice();

        // 1 thread per sequence block
        dim3 blocks( 1,1,1), threads( 1,1,1 );

        assert( offspring->getBlocksPerSequence() % 32 == 0 );

        blocks.x = offspring->getSequenceCount();
        
        if( offspring->getBlocksPerSequence() > 1024 ) {
            blocks.y = offspring->getBlocksPerSequence() / 1024 + ((offspring->getBlocksPerSequence() % 1024 ) ? 1 : 0);
            threads.x = 32;
            threads.y = 32;
        } else {
            threads.x = 32;
            threads.y = offspring->getBlocksPerSequence() / 32;
        }

        crossover_kernel<<< blocks, threads >>>( parents->getDeviceSequences(), offspring->getDeviceSequences(), sel.getDeviceList(), parents->getBlocksPerSequence(), offspring->getBlocksPerSequence() );
    }

    void computeBlockThreadDims( unsigned int M, unsigned int N, dim3 & blocks, dim3 & threads ) {
        assert( N % 32 == 0);

        unsigned int t_x = 32, t_y = N / 32;

        unsigned int b_x = M, b_y = 1;

        if( t_y > 32 ) {
            b_y = t_y / 32;

            if( t_y % 32 != 0) {
                b_y += 1;
            }

            t_y = 32;
        }

        blocks.x = b_x;
        blocks.y = b_y;

        threads.x = t_x;
        threads.y = t_y;

        std::cerr << "[ " << blocks.x << ", " << blocks.y << " ]; [ " << threads.x << ", " << threads.y << " ]" << std::endl;
    }

    void updateDevice() {
        std::cerr << "Event pool size: " << m_pool_size << " [" << m_pool_capacity << "]" << std::endl;
        std::cerr << "Event Distribution size: " << m_dist_size << " [" << m_dist_capacity << "]" << std::endl;

        assert( cudaMemcpy( m_dEventPool, m_hEventPool, m_pool_size * sizeof( event_type ), cudaMemcpyHostToDevice ) == cudaSuccess );
        assert( cudaMemcpy( m_dEventDist, m_hEventDist, m_dist_size * sizeof( unsigned int ), cudaMemcpyHostToDevice ) == cudaSuccess );
    }

    virtual ~HostCrossoverGenerator() {
        if( m_hEventDist != NULL ) {
            delete [] m_hEventDist;
            cudaFree( m_dEventDist );

            delete [] m_hEventPool;
            cudaFree( m_dEventPool );
        }
    }

protected:

    void resize( unsigned int N ) {
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

        if( 32 * N > m_pool_capacity ) {
            if( m_hEventPool != NULL ) {
                delete [] m_hEventPool;
                cudaFree( m_dEventPool );
            }

            m_hEventPool = new event_type[ 32 * N ];
            assert( cudaMalloc( (void **) &m_dEventPool, 32 * N * sizeof( event_type ) ) == cudaSuccess );

            m_pool_capacity = 32 * N;
        }
        m_pool_size = 0;
    }

    recombination_rate_parameter< double > m_recomb;

    event_type * m_hEventPool, * m_dEventPool;
    unsigned int * m_hEventDist, * m_dEventDist;

    unsigned int m_pool_size, m_pool_capacity;
    unsigned int m_dist_size, m_dist_capacity;
};

#endif  // HOST_CROSSOVER_GENERATOR_HPP_
