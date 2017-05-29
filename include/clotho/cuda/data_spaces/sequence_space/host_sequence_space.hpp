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
#ifndef HOST_SEQUENCE_SPACE_HPP_
#define HOST_SEQUENCE_SPACE_HPP_

#include "clotho/utility/bit_helper.hpp"
#include "clotho/cuda/free_space/offset_enum.hpp"
#include "clotho/utility/popcount.hpp"

template < class IntType >
class HostSequenceSpace : public clotho::utility::iStateObject {
public:
    typedef HostSequenceSpace< IntType >    self_type;
    typedef IntType                         int_type;

    typedef int_type                        block_type;

    typedef clotho::utility::BitHelper< int_type > bit_helper_type;

    HostSequenceSpace( ) :
        m_dSeqSpace(NULL)
        , m_dFreeSpace(NULL)
        , m_hFreeSpace( NULL )
        , m_blocks_per_seq(0)
        , m_allele_count(0)
        , m_seq_count(0)
        , m_capacity(0)
        , m_lost_count(0)
        , m_fixed_count(0)
        , m_free_count(0)
    {}

    void updateHost() {
        if( m_hFreeSpace != NULL ) {
            cudaMemcpy( m_hFreeSpace, m_dFreeSpace, 3 * m_blocks_per_seq * sizeof( int_type ), cudaMemcpyDeviceToHost );

            evaluateFreeSpace();
        }
    }

    void resize( unsigned int allele_count, unsigned int seq_count ) {
        unsigned int bbr = allele_count / bit_helper_type::BITS_PER_BLOCK + ((allele_count % bit_helper_type::BITS_PER_BLOCK != 0 ) ? 1 : 0);

        unsigned int new_cap = seq_count * bbr;

        if( new_cap > m_capacity ) {
            std::cerr << "Resizing Sequence Space: " << m_capacity << " -> " << new_cap << std::endl;
            if( m_dSeqSpace != NULL ) {
                cudaFree( m_dSeqSpace );
                cudaFree( m_dFreeSpace );

                delete [] m_hFreeSpace;
            }

            assert( cudaMalloc( (void **) &m_dSeqSpace, sizeof(int_type) * new_cap ) == cudaSuccess);
            assert( cudaMalloc( (void **) &m_dFreeSpace, sizeof(int_type) * bbr * 3 ) == cudaSuccess);

            m_hFreeSpace = new int_type[ 3 * bbr ];

            m_capacity = new_cap;
        }

        m_blocks_per_seq = bbr;
        m_allele_count = allele_count;
        m_seq_count = seq_count;
    }

    unsigned int getSequenceCount() const {
        return m_seq_count;
    }

    unsigned int getBlocksPerSequence() const {
        return m_blocks_per_seq;
    }

    int_type * getFixedSpace() {
        return m_hFreeSpace + FIXED_OFFSET * m_blocks_per_seq;
    }

    int_type * getLostSpace() {
        return m_hFreeSpace + LOST_OFFSET * m_blocks_per_seq;
    }

    int_type * getFreeSpace() {
        return m_hFreeSpace + FREE_OFFSET * m_blocks_per_seq;
    }

    int_type * getDeviceSequences() {
        return m_dSeqSpace;
    }

    int_type * getDeviceFreeSpace() {
        return m_dFreeSpace;
    }

    unsigned int getLostCount() const {
        return m_lost_count;
    }

    unsigned int getFreeCount() const {
        return m_free_count;
    }

    unsigned int getFixedCount() const {
        return m_fixed_count;
    }

    void get_state( boost::property_tree::ptree & state ) {
        state.put( "dimensions.rows", m_seq_count );
        state.put( "dimensions.blocks_per_row", m_blocks_per_seq);
        state.put( "dimensions.bytes_per_block", sizeof( int_type ) );
        state.put( "size", m_seq_count * m_blocks_per_seq );
        state.put( "capacity", m_capacity );
    }

    virtual ~HostSequenceSpace() {
        if( m_dSeqSpace != NULL ) {
            cudaFree( m_dSeqSpace );
            cudaFree( m_dFreeSpace );
            delete [] m_hFreeSpace;
        }
    }
protected:

    void evaluateFreeSpace() {
        m_lost_count = 0;
        for( unsigned int i = LOST_OFFSET * m_blocks_per_seq; i < (LOST_OFFSET + 1) * m_blocks_per_seq; ++i ) {
            int_type b = m_hFreeSpace[ i ];
            m_lost_count += popcount( b );
        }
        m_fixed_count = 0;
        for( unsigned int i = FIXED_OFFSET * m_blocks_per_seq; i < (FIXED_OFFSET + 1) * m_blocks_per_seq; ++i ) {
            int_type b = m_hFreeSpace[ i ];
            m_fixed_count += popcount( b );
        }

        m_free_count = 0;
        for( unsigned int i = FREE_OFFSET * m_blocks_per_seq; i < (FREE_OFFSET + 1) * m_blocks_per_seq; ++i ) {
            int_type b = m_hFreeSpace[ i ];
            m_free_count += popcount( b );
        }
    }

    int_type    * m_dSeqSpace;

    int_type    * m_hFreeSpace, * m_dFreeSpace;

    unsigned int m_blocks_per_seq, m_seq_count, m_allele_count;
    size_t      m_capacity;

    unsigned int m_lost_count, m_fixed_count, m_free_count;
};
#endif  // HOST_SEQUENCE_SPACE_HPP_
