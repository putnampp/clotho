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
#ifndef HOST_ALLELE_SPACE_HPP_
#define HOST_ALLELE_SPACE_HPP_

template < class RealType >
class HostAlleleSpace : public clotho::utility::iStateObject {
public:
    typedef HostAlleleSpace< RealType > self_type;

    typedef RealType                    real_type;

    typedef real_type                   location_type;
    typedef unsigned int                age_type;

    HostAlleleSpace() :
        m_dLoc( NULL )
        , m_hLoc( NULL )
        , m_size(0)
        , m_capacity(0)
        , m_dSize(0)
        , m_dCapacity(0)
    {}

    void updateDevice() {
        if( m_capacity > m_dCapacity ) {
            if( m_dLoc != NULL ) {
                std::cerr << "Free device locations" << std::endl;
                assert( cudaFree( m_dLoc ) == cudaSuccess);
            }
            assert( cudaMalloc((void **)&m_dLoc, m_capacity * sizeof(location_type)) == cudaSuccess );
            m_dCapacity = m_capacity;
        }

        m_dSize = m_size;
        assert( m_dSize <= m_dCapacity );
        assert( cudaMemcpy( m_dLoc, m_hLoc, sizeof( location_type ) * m_dCapacity, cudaMemcpyHostToDevice ) == cudaSuccess);
    }

    void get_state( boost::property_tree::ptree & state ) {
        boost::property_tree::ptree l, a;

        for( unsigned int i = 0; i < m_size; ++i ) {
            clotho::utility::add_value_array( l, m_hLoc[ i ] );
            clotho::utility::add_value_array( a, m_hAge[ i ] );
        }

        state.put( "size", m_size );
        state.put( "capacity", m_capacity );
        state.put_child( "locations", l);
        state.put_child( "age", a);
    }

    size_t getDeviceAlleleCount() const {
        return m_dSize;
    }

    size_t getDeviceMaxAlleles() const {
        return m_dCapacity;
    }

    size_t getAlleleCount() const {
        return m_size;
    }

    size_t getMaxAlleleCount() const {
        return m_capacity;
    }

    void update( unsigned int idx, location_type l, age_type a ) {
        assert( idx < m_size );

        this->m_hLoc[ idx ] = l;
        this->m_hAge[ idx ] = a;
    }

    void update( unsigned int a_idx, self_type & others, unsigned int idx ) {
        assert( a_idx < m_size );
        assert( idx < others.m_size );

        this->m_hLoc[ a_idx ] = others.m_hLoc[ idx ];
        this->m_hAge[ a_idx ] = others.m_hAge[ idx ];
    }

    void resize( unsigned int N ) {
        size_t new_cap = N % 1024;
        if( new_cap > 0 ) {
            new_cap = 1024 - new_cap;
        }

        new_cap += N + 1024;

        if( new_cap > m_capacity ) {
            std::cerr << "Resizing alleles: " << m_capacity << " -> " << new_cap << std::endl;
            location_type * tmp_loc = new location_type[ new_cap ];
            age_type * tmp_age = new age_type[ new_cap ];
            
            if( m_hLoc != NULL ) {
                memcpy( tmp_loc, m_hLoc, sizeof( location_type ) * m_size);
                memcpy( tmp_age, m_hAge, sizeof( age_type ) * m_size);

                delete [] m_hLoc;
                delete [] m_hAge;
            }

            m_hLoc = tmp_loc;
            m_hAge = tmp_age;

            m_capacity = new_cap;
        }

        m_size = N;
    }

    location_type * getDeviceLocations() {
        return m_dLoc;
    }

    virtual ~HostAlleleSpace() {
        if( m_hLoc != NULL ) {
            delete [] m_hLoc;
            delete [] m_hAge;
        }

        if( m_dLoc != NULL ) {
            cudaFree( m_dLoc );
        }
    }

protected:

    location_type * m_hLoc, * m_dLoc;

    age_type * m_hAge;

    size_t m_size, m_capacity;
    size_t m_dSize, m_dCapacity;
};

#endif  // HOST_ALLELE_SPACE_HPP_
