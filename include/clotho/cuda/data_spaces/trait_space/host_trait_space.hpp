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
#ifndef HOST_TRAIT_SPACE_HPP_
#define HOST_TRAIT_SPACE_HPP_

#include "clotho/data_spaces/phenotype_evaluator/trait_count_parameter.hpp"

template < class RealType >
class HostTraitSpace : public clotho::utility::iStateObject {
public:

    typedef HostTraitSpace< RealType >  self_type;

    typedef RealType                    real_type;
    typedef real_type                   weight_type;

    HostTraitSpace( boost::property_tree::ptree & config ) :
        m_hTraitSpace( NULL )
        , m_dTraitSpace( NULL )
        , m_allele_count(0)
        , m_trait_count(0)
        , m_capacity(0)
        , m_dCapacity(0)
    {
        trait_count_parameter param( config );
        m_trait_count = param.m_trait_count;
    }

    void get_state( boost::property_tree::ptree & s ) {
        boost::property_tree::ptree d;
        for( unsigned int i = 0; i < m_trait_count; ++i ) {
            boost::property_tree::ptree t;
            for( unsigned int j = 0; j < m_allele_count; ++j ) {
                clotho::utility::add_value_array( t, m_hTraitSpace[ i * m_allele_count + j ] );
            }
            d.push_back( std::make_pair( "", t ));
        }

        s.put( "trait_count", m_trait_count);
        s.put( "allele_count", m_allele_count );
        s.put( "size", m_trait_count * m_allele_count );
        s.put( "capacity", m_capacity );
        s.put_child( "data", d );
    }


    weight_type * getDeviceWeights() {
        return m_dTraitSpace;
    }

    size_t getTraitCount() const {
        return m_trait_count;
    }

    size_t getAlleleCount() const {
        return m_allele_count;
    }

    void update( unsigned int allele_idx, unsigned int trait_idx, weight_type w ) {
        assert( allele_idx < m_allele_count );
        assert( trait_idx < m_trait_count );

        m_hTraitSpace[ trait_idx * m_allele_count + allele_idx ] = w;
    }

    void updateDevice() {
        resizeDevice();
        assert( cudaMemcpy( m_dTraitSpace, m_hTraitSpace, m_dCapacity * sizeof(weight_type), cudaMemcpyHostToDevice ) == cudaSuccess );
    }

    void resizeDevice() {
        if( m_capacity > m_dCapacity ) {
            if( m_dTraitSpace != NULL ) {
                assert( cudaFree( m_dTraitSpace ) == cudaSuccess );
            }
            assert( cudaMalloc((void **) &m_dTraitSpace, m_capacity * sizeof(weight_type)) == cudaSuccess);
            m_dCapacity = m_capacity;
        }
    }

    void updateDevice( cudaStream_t & stream ) {
        resizeDevice();
        assert( cudaMemcpyAsync( m_dTraitSpace, m_hTraitSpace, m_dCapacity * sizeof(weight_type), cudaMemcpyHostToDevice, stream ) == cudaSuccess );
    }

    void update( unsigned int a_idx, self_type & other, unsigned int b_idx ) {
        assert( a_idx < m_allele_count );
        assert( b_idx < other.m_allele_count );
        assert( m_trait_count <= other.m_trait_count );

        for( unsigned int i = 0; i < m_trait_count; ++i ) {
            m_hTraitSpace[ a_idx ] = other.m_hTraitSpace[ b_idx ];
            a_idx += m_allele_count;
            b_idx += other.m_allele_count;
        }
    }

    void resize( HostAlleleSpace < RealType > & alleles ) {
        resize( alleles.getMaxAlleleCount(), m_trait_count );
    }

    virtual ~HostTraitSpace() {
        if( m_dTraitSpace != NULL ) {
            delete [] m_hTraitSpace;
        }
        if( m_dTraitSpace != NULL ) {
            cudaFree( m_dTraitSpace );
        }
    }

protected:

    void resize( unsigned int allele_count, unsigned int trait_count ) {
        size_t new_cap = allele_count * trait_count;
        if( new_cap > m_capacity ) {

            weight_type * tmp_weight = new weight_type[ new_cap ];
            
            if( m_dTraitSpace != NULL )  {

                // copy existing weights over
                for( unsigned int i = 0; i < m_trait_count; ++i ) {
                    memcpy( tmp_weight + i * allele_count, m_hTraitSpace + i * m_allele_count, sizeof( weight_type ) * m_allele_count );
                }

                delete [] m_hTraitSpace;
            }

            m_hTraitSpace = tmp_weight;

            m_capacity = new_cap;
        }

        m_allele_count = allele_count;
        m_trait_count = trait_count;
    }

    weight_type   * m_hTraitSpace;
    weight_type   * m_dTraitSpace;

    size_t    m_allele_count, m_trait_count;

    size_t    m_capacity, m_dCapacity;
};

#endif  // HOST_TRAIT_SPACE_HPP_
