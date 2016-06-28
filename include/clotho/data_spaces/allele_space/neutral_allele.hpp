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
#ifndef CLOTHO_NEUTRAL_ALLELE_HPP_
#define CLOTHO_NEUTRAL_ALLELE_HPP_

#include "clotho/data_spaces/allele_space/base_allele.hpp"

#include "clotho/utility/state_object.hpp"
#include "clotho/utility/log_helper.hpp"

#include "clotho/utility/bit_helper.hpp"

namespace clotho {
namespace genetics {

template < class PositionType >
class neutral_allele_vectorized : public base_allele_vectorized< PositionType > {
public:
    typedef base_allele_vectorized< PositionType >      base_type;
    typedef neutral_allele_vectorized< PositionType >   self_type;

    typedef unsigned long long                          block_type;
    typedef block_type *                                neutrality_vector_type;

    typedef clotho::utility::BitHelper< block_type >    bit_helper_type;

    neutral_allele_vectorized( size_t a = 0 ) :
        base_type( a )
        , m_neutral( NULL )
        , m_block_size(0)
        , m_block_alloc(0)
    {
        this->grow( a );
    }

    bool getNeutralAt( size_t index ) const {
        assert( 0 <= index && index < this->m_pos_size );

        size_t bidx = index / bit_helper_type::BITS_PER_BLOCK;
        block_type bmask = bit_helper_type::bit_offset(index);

        return m_neutral[ bidx ] & bmask;
    }

    void setNeutralAt( size_t index, bool neu ) {
        assert( 0 <= index && index < this->m_pos_size );

        size_t bidx = index / bit_helper_type::BITS_PER_BLOCK;
        block_type mask = bit_helper_type::bit_offset(index);

        if( ((m_neutral[ bidx ] & mask) == mask) != neu )
            m_neutral[ bidx ] ^= mask;
    }

//    bool getNeutralAt( size_t index ) const {
//        return m_neutral[ index ];
//    }
//
//    void setNeutralAt( size_t index, bool neu ) {
//        if( index >= base_type::size() ) {
//            resize( index + 1 );
//        }
//
//        m_neutral[ index ] = neu;
//    }

//    bool isAllNeutral() const {
//        const_neutrality_iterator first = m_neutral.begin(), last = this->m_neutral.end();
//        bool all_neutral = true;
//        while( all_neutral && first != last ) {
//            all_neutral = *first++;
//        }
//
//        return all_neutral;
//    }
    bool isAllNeutral() const {
        for( size_t i = 0; i < m_block_size - 1; ++i ) {
            if( m_neutral[i] != bit_helper_type::ALL_SET )
                return false;
        }

        block_type mask = bit_helper_type::low_bit_mask( this->m_pos_size - 1 );

//        std::cerr << "Comparing: " << std::hex << m_neutral[m_block_size - 1] << " to " << mask << std::dec << std::endl;

//        std::cerr << "Checking: " << (m_block_size - 1) << std::endl;

        return m_neutral[m_block_size - 1] == mask;
    }

    void inherit( self_type & parent ) {
        // inherit positions
        assert( parent.size() <= this->size() );
        memcpy( this->m_positions, parent.m_positions, parent.size() * sizeof(typename base_type::base_type::position_type) );

        // inherit neutral
        assert( parent.m_block_size < this->m_block_alloc );
        std::copy( parent.m_neutral, parent.m_neutral + parent.m_block_size, this->m_neutral );
    }

//    neutrality_iterator neutral_begin() {
//        return m_neutral.begin();
//    }
//
//    neutrality_iterator neutral_end() {
//        return m_neutral.end();
//    }
//
//    const_neutrality_iterator neutral_begin() const {
//        return m_neutral.begin();
//    }
//
//    const_neutrality_iterator neutral_end() const {
//        return m_neutral.end();
//    }

    void push_back( self_type & other, size_t idx ) {
        size_t e = this->size();

        this->resize( e + 1 );

        this->setPositionAt( e, other.getPositionAt( idx ) );
        this->setNeutralAt( e, other.getNeutralAt( idx ) );
    }

    virtual size_t grow( size_t rows ) {
        this->resize( rows );

        return this->size();
    }

    virtual ~neutral_allele_vectorized() {
        if( m_neutral != NULL ) {
            delete [] m_neutral;
        }
    }

protected:

    virtual void resize( size_t s ) {
        base_type::resize( s );

        size_t bc = bit_helper_type::padded_block_count( s );
        if( bc > m_block_alloc ) {
            if( m_neutral != NULL ) {
                delete [] m_neutral;
            }

            m_neutral = new block_type[ bc ];
            m_block_alloc = bc;

        }

        memset( m_neutral, 0, m_block_alloc * sizeof(block_type));
        m_block_size = bc;
    }

    neutrality_vector_type  m_neutral;
    size_t m_block_size, m_block_alloc;
};

}   // namespace genetics
}   // namespace clotho

namespace clotho {
namespace utility {

template < class PositionType >
struct state_getter< clotho::genetics::neutral_allele_vectorized< PositionType > > {
    typedef clotho::genetics::neutral_allele_vectorized< PositionType > object_type;

    void operator()( boost::property_tree::ptree & s, object_type & obj ) {
        size_t all_count = obj.allele_count();
        size_t i = 0;
        while( i < all_count ) {
            boost::property_tree::ptree all;
            all.put( "position", obj.getPositionAt(i) );
            all.put( "neutral", obj.getNeutralAt(i) );

            s.push_back( std::make_pair("", all ) );
            ++i;
        }
        
    }
};

}   // namespace utility
}   // namespace clotho

#endif  // CLOTHO_NEUTRAL_ALLELE_HPP_
