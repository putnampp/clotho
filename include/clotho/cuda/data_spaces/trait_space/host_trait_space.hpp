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
    {}

    void get_state( boost::property_tree::ptree & s ) {

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
        assert( cudaMemcpy( m_dTraitSpace, m_hTraitSpace, m_allele_count * m_trait_count * sizeof( weight_type ), cudaMemcpyHostToDevice ) == cudaSuccess );
    }

    void resize( HostAlleleSpace < RealType > & alleles ) {
        resize( alleles.getMaxAlleleCount(), m_trait_count );
    }

    virtual ~HostTraitSpace() {
        if( m_dTraitSpace != NULL ) {
            cudaFree( m_dTraitSpace );
            delete [] m_hTraitSpace;
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

                cudaFree( m_dTraitSpace );

                delete [] m_hTraitSpace;
            }

            m_hTraitSpace = tmp_weight;

            assert( cudaMalloc( (void **) &m_dTraitSpace, sizeof( real_type ) * new_cap ) == cudaSuccess);
            m_capacity = new_cap;
        }

        m_allele_count = allele_count;
        m_trait_count = trait_count;
    }

    real_type   * m_hTraitSpace;
    real_type   * m_dTraitSpace;

    size_t    m_allele_count, m_trait_count;

    size_t    m_capacity;
};

#endif  // HOST_TRAIT_SPACE_HPP_
