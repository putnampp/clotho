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
#ifndef HOST_PHENOTYPE_SPACE_HPP_
#define HOST_PHENOTYPE_SPACE_HPP_

template < class RealType >
class HostPhenotypeSpace : public clotho::utility::iStateObject {
public:
    typedef HostPhenotypeSpace< RealType >  self_type;

    typedef RealType                        real_type;

    template < class IntType >
    void operator()( HostSequenceSpace< IntType > & seq_space, HostTraitSpace< RealType > & trait_space ) {
        resize( seq_space.getSequenceCount(), trait_space.getTraitCount() );
    }

    unsigned int getSequenceCount() const {
        return m_seq_count;
    }

    unsigned int getTraitCount() const {
        return m_trait_count;
    }

    real_type * getDevicePhenotype() {
        return m_dPhenoSpace;
    }

    virtual HostPhenotypeSpace() {
        if( m_dPhenoSpace != NULL ) {
            cudaFree( m_dPhenoSpace );
        }
    }
 
protected:
    void resize( unsigned int seq_count, unsigned int trait_count ) {
        size_t new_cap = seq_count * trait_count;

        if( new_cap > m_capacity ) {
            if( m_dPhenoSpace != NULL ) {
                cudaFree( m_dPhenoSpace );
            }

            cudaMalloc( (void **) m_dPhenoSpace, sizeof( real_type ) * new_cap );

            m_capacity = new_cap;
        }

        m_seq_count = seq_count;
        m_trait_count = trait_count;
    }

    real_type   * m_dPhenoSpace;

    unsigned int m_seq_count, m_trait_count;
    size_t       m_capacity;
};

#endif  // HOST_PHENOTYPE_SPACE_HPP_
