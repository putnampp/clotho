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
#ifndef HOST_MUTATE_GENERATOR_HPP_
#define HOST_MUTATE_GENERATOR_HPP_

#include <boost/random/poisson_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/normal_distribution.hpp>

#include "clotho/mutation/mutation_rate_parameter.hpp"
#include "clotho/data_spaces/generators/weight_parameter.hpp"
#include "clotho/data_spaces/generators/neutral_parameter.hpp"

#include "clotho/cuda/mutation/mutate_kernel.hpp"

class HostMutateGenerator {
public:

    HostMutateGenerator( boost::property_tree::ptree & config ) :
        m_rate( config )
        , m_neutral(config)
        , m_hAlleleIndexPool(NULL)
        , m_dAlleleIndexPool(NULL)
        , m_hSeqDist(NULL)
        , m_dSeqDist(NULL)
        , m_pool_size(0)
        , m_pool_capacity(0)
        , m_dist_size(0)
        , m_dist_capacity(0)
    {}

    template < class RNG, class RealType, class IntType >
    unsigned int initialize( RNG & eng, HostPopulationSpace< RealType, IntType > * parents, unsigned int N, unsigned int all_count ) {
        double lambda = N * m_rate.m_mu;
        boost::random::poisson_distribution< unsigned int, double > dist( lambda );

        unsigned int M = dist( eng );

        std::cerr << "New mutations: " << M << std::endl;

        resize( N, M );

        return buildAllele( parents, all_count);
    }

    template < class RNG, class RealType >
    void generate( RNG & eng, HostAlleleSpace< RealType > & alleles, HostTraitSpace< RealType > & traits, unsigned int age ) {

        updateSpace( eng, alleles, traits, age );
    }

    template < class RealType, class IntType >
    void operator()( HostPopulationSpace< RealType, IntType > * pop ) {
        updateDevice();

        dim3 blocks( pop->getSequenceCount(), 1, 1 ), threads( 1,1,1 );

        std::cerr << "Mutate: [" << blocks.x << ", " << blocks.y << "]; [" << threads.x << ", " << threads.y << "]" << std::endl;
        assert( pop->getDeviceSequences() != NULL );
        assert( m_dAlleleIndexPool != NULL );
        assert( m_dSeqDist != NULL );
        mutate<<< blocks, threads >>>( pop->getDeviceSequences(), m_dAlleleIndexPool, m_dSeqDist, pop->getBlocksPerSequence() );
    }

    void updateDevice() {
        assert( cudaMemcpy( m_dSeqDist, m_hSeqDist, sizeof( unsigned int ) * m_dist_size, cudaMemcpyHostToDevice ) == cudaSuccess );

        std::cerr << "Sequence Distribution: ";
        for( unsigned int i = 0; i < 10; ++i ) {
            std::cerr << m_hSeqDist[ i ] << ", ";
        }
        std::cerr << " ... , " << m_hSeqDist[ m_dist_size - 1] << std::endl;
        assert( cudaMemcpy( m_dAlleleIndexPool, m_hAlleleIndexPool, sizeof( unsigned int ) * m_pool_size, cudaMemcpyHostToDevice ) == cudaSuccess );

        std::cerr << "Allele Pool: ";
        for( unsigned int i = 0; i < 10; ++i ) {
            std::cerr << m_hAlleleIndexPool[ i ] << ", ";
        }
        std::cerr << " ... , " << m_hAlleleIndexPool[ m_pool_size - 1] << std::endl;
    }

    virtual ~HostMutateGenerator() {
        if( m_hAlleleIndexPool != NULL ) {
            delete [] m_hAlleleIndexPool;
            delete [] m_hSeqDist;

            cudaFree( m_dAlleleIndexPool );
            cudaFree( m_dSeqDist );
        }
    }

protected:

    template < class RealType, class IntType >
    unsigned int buildAllele( HostPopulationSpace< RealType, IntType > * parents, unsigned int all_count ) {
        typedef HostPopulationSpace< RealType, IntType > population_space_type;
        typedef typename population_space_type::block_type    block_type;
        typedef typename population_space_type::sequence_space_type::bit_helper_type bit_helper_type;

        unsigned int mut_idx = 0;
        for( unsigned int i = 0; i < parents->getBlocksPerSequence() && mut_idx < m_pool_size; ++i ) {
            block_type b = parents->getFreeSpace()[ i ];

            unsigned int j = i * bit_helper_type::BITS_PER_BLOCK;
            while( b > 0 ) {
                if( b & 1) {
                    m_hAlleleIndexPool[ mut_idx++ ] = j;
                    if( mut_idx >= m_pool_size ) break;
                }
                b >>= 1;
                ++j;
            }
        }

        while( mut_idx < m_pool_size ) {
            m_hAlleleIndexPool[ mut_idx++ ] = all_count++;
        }

        return all_count;
    }

    template < class RNG, class RealType >
    void updateSpace( RNG & eng, HostAlleleSpace< RealType > & alleles, HostTraitSpace< RealType > & traits, unsigned int age ) {
        typedef HostAlleleSpace< RealType > allele_space_type;
        typedef typename allele_space_type::location_type location_type;

        typedef HostTraitSpace< RealType > trait_space_type;
        typedef typename trait_space_type::weight_type weight_type;

        typedef boost::random::uniform_01< location_type > location_distribution_type;
        typedef boost::random::uniform_int_distribution< unsigned int > sequence_distribution_type;
        typedef boost::random::normal_distribution< weight_type >       weight_distribution_type;

        location_distribution_type loc_dist;
        sequence_distribution_type seq_dist( 0, m_dist_size - 2);
        weight_distribution_type weight_dist( 0.0, 1.0);

        for( unsigned int i = 0; i < m_pool_size; ++i ) {
            unsigned int idx = m_hAlleleIndexPool[ i ];
            location_type loc = loc_dist( eng );

            if( idx >= alleles.getAlleleCount() ) {
                alleles.push_back( loc, age );
            } else {
                alleles.update( idx, loc, age );
            }

            for( unsigned int t = 0; t < traits.getTraitCount(); ++t ) {
                weight_type w = 0.0;

                if( !m_neutral.m_dist( eng ) ) {
                    w = weight_dist( eng );
                }

                traits.update( idx, t, w );
            }

             m_hSeqDist[ seq_dist( eng ) ] += 1;
        }

        // scan right of sequence distribution
        unsigned int N = m_hSeqDist[ 0 ];
        m_hSeqDist[ 0 ] = 0;
        for( unsigned int i = 1; i < m_dist_size; ++i ) {
            unsigned int tmp = m_hSeqDist[ i ];
            m_hSeqDist[ i ] = N;
            N += tmp;
        }
    }

    void resize( unsigned int seq_count, unsigned int mut_count ) {
        if( mut_count > m_pool_capacity ) {
            if( m_hAlleleIndexPool != NULL ) {
                delete [] m_hAlleleIndexPool;
                cudaFree( m_dAlleleIndexPool );
            }

            m_hAlleleIndexPool = new unsigned int[ mut_count ];

            assert( cudaMalloc( (void **) &m_dAlleleIndexPool, sizeof(unsigned int) * mut_count ) == cudaSuccess );

            m_pool_capacity = mut_count;
        }

        m_pool_size = mut_count;

        seq_count += 1;
        if( seq_count > m_dist_capacity ) {
            if( m_hSeqDist != NULL ) {
                delete [] m_hSeqDist;
                cudaFree( m_dSeqDist );
            }

            m_hSeqDist = new unsigned int[ seq_count ];
            assert( cudaMalloc( (void **) & m_dSeqDist, sizeof(unsigned int) * seq_count ) == cudaSuccess );

            m_dist_capacity = seq_count;
        }

        m_dist_size = seq_count;

        memset( m_hSeqDist, 0, sizeof( unsigned int ) * m_dist_size );
    }

    mutation_rate_parameter< double > m_rate;
    clotho::genetics::neutral_parameter< double > m_neutral;

    unsigned int * m_hAlleleIndexPool, * m_dAlleleIndexPool;
    unsigned int * m_hSeqDist, * m_dSeqDist;

    size_t m_pool_size, m_pool_capacity;
    size_t m_dist_size, m_dist_capacity;
};

#endif  // HOST_MUTATE_GENERATOR_HPP_
