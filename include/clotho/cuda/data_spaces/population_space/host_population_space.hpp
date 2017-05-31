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
#ifndef HOST_POPULATION_SPACE_DEF_HPP_
#define HOST_POPULATION_SPACE_DEF_HPP_

#include "clotho/cuda/data_spaces/sequence_space/host_sequence_space.hpp"
#include "clotho/cuda/data_spaces/phenotype_space/host_phenotype_space.hpp"
#include "clotho/cuda/data_spaces/allele_space/host_allele_space.hpp"

template < class RealType, class IntType >
class HostPopulationSpace : public clotho::utility::iStateObject {
public:
    typedef HostPopulationSpace< RealType, IntType > self_type;

    typedef HostSequenceSpace< IntType >                sequence_space_type;
    typedef HostPhenotypeSpace< RealType >              phenotype_space_type;

    typedef typename sequence_space_type::block_type    block_type;

    typedef typename phenotype_space_type::phenotype_type  phenotype_type;

    typedef real_type                                   fitness_type;

    HostPopulationSpace( ) :
        m_dFitness( NULL )
        , m_hFitness( NULL )
        , m_size(0)
        , m_capacity(0)
    {}

    fitness_type * getDeviceFitness() {
        return m_dFitness;
    }

    fitness_type * getHostFitness() {
        return m_hFitness;
    }

    void updateHost() {
        m_seqs.updateHost();
        assert( cudaMemcpy( m_hFitness, m_dFitness, m_capacity * sizeof( fitness_type ), cudaMemcpyDeviceToHost ) == cudaSuccess);
    }

    void updateHostFitnessAsync( cudaStream_t & stream ) {
        assert( cudaMemcpyAsync( m_hFitness, m_dFitness, m_capacity * sizeof( fitness_type ), cudaMemcpyDeviceToHost, stream ) == cudaSuccess);
    }

    void evaluateFreeSpace() {
        m_seqs.evaluateFreeSpace();
    }

    void updateHostFreeSpaceAsync( cudaStream_t & stream) {
        m_seqs.updateHost( stream );
    }

    block_type * getLostSpace() {
        return m_seqs.getLostSpace();
    }

    block_type * getFixedSpace() {
        return m_seqs.getFixedSpace();
    }

    block_type * getFreeSpace() {
        return m_seqs.getFreeSpace();
    }

    unsigned int getIndividualCount() const {
        return m_size;
    }

    void resize( HostAlleleSpace< RealType > & allele_space , HostTraitSpace< RealType > & trait_space, unsigned int nSeqs ) {
        m_seqs.resize( allele_space.getMaxAlleleCount(), nSeqs );
        m_phenos.resize( nSeqs, trait_space.getTraitCount() );

        unsigned int fit_count = m_phenos.getSequenceCount() / 2;

        assert( m_phenos.getSequenceCount() % 2 == 0 );

        if( fit_count > m_capacity ) {
            if( m_dFitness != NULL ) {
                assert( cudaFree( m_dFitness ) == cudaSuccess);
                delete [] m_hFitness;
            }

            assert( cudaMalloc( (void **) &m_dFitness, sizeof( fitness_type ) * fit_count ) == cudaSuccess );

            m_hFitness = new fitness_type[ fit_count ];

            m_capacity = fit_count; 
        }

        m_size = fit_count;
    }

    block_type * getDeviceSequences() {
        return m_seqs.getDeviceSequences();
    }

    unsigned int getBlocksPerSequence() const {
        return m_seqs.getBlocksPerSequence();
    }

    unsigned int getSequenceCount() const {
        return m_seqs.getSequenceCount();
    }

    unsigned int getTraitCount() const {
        return m_phenos.getTraitCount();
    }

    phenotype_type * getDevicePhenotypes() {
        return m_phenos.getDevicePhenotypes();
    }

    block_type * getDeviceFreeSpace() {
        return m_seqs.getDeviceFreeSpace();
    }

    unsigned int getFreeCount() const {
        return m_seqs.getFreeCount();
    }

    unsigned int getFixedCount() const {
        return m_seqs.getFixedCount();
    }

    unsigned int getLostCount() const {
        return m_seqs.getLostCount();
    }

    void get_state( boost::property_tree::ptree & state ) {
        updateHost();

        boost::property_tree::ptree s, p, f, d;
        m_seqs.get_state( s );
        m_phenos.get_state( p );

        for( unsigned int i = 0; i < m_size; ++i ) {
            clotho::utility::add_value_array( d, m_hFitness[ i ] );
        }

        f.put( "size", m_size );
        f.put( "capacity", m_capacity );
        f.put_child( "data", d );

        state.put_child( "sequence_space", s );
        state.put_child( "phenotype_space", p );
        state.put_child( "fitness", f );
    }

    virtual ~HostPopulationSpace() { 
        if( m_dFitness != NULL ) {
            cudaFree( m_dFitness );
            delete [] m_hFitness;
        }
    }

protected:

    sequence_space_type m_seqs;
    phenotype_space_type m_phenos;

    fitness_type * m_dFitness, * m_hFitness;

    size_t  m_size, m_capacity;
};

#endif  // HOST_POPULATION_SPACE_DEF_HPP_

