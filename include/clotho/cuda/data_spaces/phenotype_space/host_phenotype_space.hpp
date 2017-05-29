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
    typedef real_type                       phenotype_type;

    HostPhenotypeSpace() :
        m_hPhenoSpace( NULL )
        , m_dPhenoSpace( NULL )
        , m_seq_count(0)
        , m_trait_count(0)
        , m_capacity(0)
    {}

    unsigned int getSequenceCount() const {
        return m_seq_count;
    }

    unsigned int getTraitCount() const {
        return m_trait_count;
    }

    phenotype_type * getDevicePhenotypes() {
        return m_dPhenoSpace;
    }

    void get_state( boost::property_tree::ptree & s ) {
        updateHost();

        boost::property_tree::ptree d;
        for( unsigned int i = 0; i < m_trait_count; ++i ) {
            boost::property_tree::ptree t;
            for( unsigned int j = 0; j < m_seq_count; ++j ) {
                clotho::utility::add_value_array( t, m_hPhenoSpace[ i * m_trait_count + j ] );
            }
            d.push_back( std::make_pair( "", t ) );
        }
        
        s.put( "size", m_seq_count * m_trait_count );
        s.put( "capacity", m_capacity );
        s.put_child( "data", d );
    }

    void updateHost() {
        assert( cudaMemcpy( m_hPhenoSpace, m_dPhenoSpace, m_seq_count * m_trait_count * sizeof( phenotype_type ), cudaMemcpyDeviceToHost ) == cudaSuccess );
    }

    virtual ~HostPhenotypeSpace() {
        if( m_dPhenoSpace != NULL ) {
            cudaFree( m_dPhenoSpace );
        }
    }
 
    void resize( unsigned int seq_count, unsigned int trait_count ) {
        size_t new_cap = seq_count * trait_count;

        if( new_cap > m_capacity ) {
            if( m_hPhenoSpace != NULL ) {
                delete [] m_hPhenoSpace;
                cudaFree( m_dPhenoSpace );
            }

            assert( cudaMalloc( (void **) &m_dPhenoSpace, sizeof( phenotype_type ) * new_cap ) == cudaSuccess );

            m_hPhenoSpace = new phenotype_type[ new_cap ];
            m_capacity = new_cap;
        }

        m_seq_count = seq_count;
        m_trait_count = trait_count;
    }

protected:
    real_type   *m_hPhenoSpace, * m_dPhenoSpace;

    unsigned int m_seq_count, m_trait_count;
    size_t       m_capacity;
};

#endif  // HOST_PHENOTYPE_SPACE_HPP_
