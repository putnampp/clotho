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
#ifndef CLOTHO_ALLELE_SPACE_VECTOR_HPP_
#define CLOTHO_ALLELE_SPACE_VECTOR_HPP_

#include <vector>
#include <boost/dynamic_bitset.hpp>

#include "clotho/utility/state_object.hpp"
#include "clotho/utility/log_helper.hpp"

#include "clotho/utility/bit_helper.hpp"

namespace clotho {
namespace genetics {

template < class PositionType, class SizeType = unsigned int >
class AlleleSpace {
public:

    typedef PositionType                    position_type;
    typedef std::vector< position_type >    position_vector;

    typedef SizeType                        size_type;

    typedef boost::dynamic_bitset< unsigned int >   neutral_vector;

    typedef unsigned int                age_type;
    typedef std::vector< age_type >     age_vector;

    typedef AlleleSpace< PositionType, SizeType >   self_type;

    static const unsigned int ALLELES_PER_BLOCK = clotho::utility::BitHelper< unsigned int >::BITS_PER_BLOCK; 

    AlleleSpace() {}

    position_vector & getPositions() {
        return m_positions;
    }

    position_type getAllelePosition( size_type index ) const {
        return m_positions[index];
    }

    bool getAlleleNeutral( size_type index ) const {
        return m_neutral[ index ];
    }

    age_type getAlleleAge( size_type index ) const {
        return m_age[ index ];
    }

    bool isAllNeutral() const {
        return m_neutral.all();
    }

    void setAllele( size_type index, position_type pos, bool is_neutral, unsigned int age ) {
#ifdef DEBUGGING
        assert( index < m_positions.size() );
#endif // DEBUGGING

        m_positions[ index ] = pos;
        m_neutral.set( index, is_neutral );
        m_age[ index ] = age;
    }

//    void setAllele( size_type index, position_type pos, bool is_neutral, unsigned int age ) {
//        if( index < m_positions.size() ) {
//            m_positions[ index ] = pos;
//            m_neutral.set( index, is_neutral );
//            m_age[ index ] = age;
//        } else {
//            do {
//                m_positions.push_back( pos );
//                m_neutral.push_back( is_neutral );
//                m_age.push_back(age);
//            } while( m_positions.size() <= index );
//#ifdef DEBUGGING
//            assert( m_positions.size() == m_neutral.size() );
//            assert( m_positions.size() == m_age.size() );
//#endif  // DEBUGGING
//        }
//    }

    void append( const self_type & other, unsigned int idx ) {
#ifdef DEBUGGING
        assert( idx < other.size() );
#endif  // DEBUGGING

        m_positions.push_back( other.m_positions[ idx ] );
        m_neutral.push_back( other.m_neutral[ idx ] );
        m_age.push_back( other.m_age[ idx ] );
    }

    size_t size() const {
        return m_positions.size();
    }

    void grow() {
        m_positions.push_back( 0 );
        m_neutral.push_back( false );
        m_age.push_back( 0 );
#ifdef DEBUGGING
        assert( m_positions.size() == m_neutral.size() );
        assert( m_positions.size() == m_age.size() );
#endif  // DEBUGGING
    }

    virtual ~AlleleSpace() {}
protected:

    position_vector m_positions;
    neutral_vector  m_neutral;
    age_vector      m_age;
};

}   // namespace genetics
}   // namespace clotho


namespace clotho {
namespace utility {

template < class PositionType, class SizeType >
struct state_getter< clotho::genetics::AlleleSpace< PositionType, SizeType > > {
    typedef clotho::genetics::AlleleSpace< PositionType, SizeType > object_type;

    void operator()( boost::property_tree::ptree & s, object_type & obj ) {
        size_t all_count = obj.size();
        size_t i = 0;
        while( i < all_count ) {
            boost::property_tree::ptree all;
            all.put( "position", obj.getAllelePosition(i) );
            all.put( "neutral", obj.getAlleleNeutral(i) );
            all.put( "age", obj.getAlleleAge(i) );

            s.push_back( std::make_pair("", all ) );
            ++i;
        }

    }
};

}   // namespace utility
}   // namespace clotho
#endif  // CLOTHO_ALLELE_SPACE_VECTOR_HPP_
