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

namespace clotho {
namespace genetics {

template < class PositionType >
class neutral_allele_vectorized : public base_allele_vectorized< PositionType > {
public:
    typedef base_allele_vectorized< PositionType >      base_type;
    typedef neutral_allele_vectorized< PositionType >   self_type;

    typedef std::vector< bool >                             neutrality_vector_type;
    typedef typename neutrality_vector_type::iterator       neutrality_iterator;
    typedef typename neutrality_vector_type::const_iterator const_neutrality_iterator;

    neutral_allele_vectorized( size_t a = 0 ) :
        base_type( a )
    {
        this->grow( a );
    }

    bool getNeutralAt( size_t index ) const {
        return m_neutral[ index ];
    }

    void setNeutralAt( size_t index, bool neu ) {
        if( index >= base_type::size() ) {
            resize( index + 1 );
        }

        m_neutral[ index ] = neu;
    }

    void inherit( self_type & parent ) {
        std::copy( parent.position_begin(), parent.position_end(), this->position_begin() );
        std::copy( parent.neutral_begin(), parent.neutral_end(), this->m_neutral.begin() );
    }

    neutrality_iterator neutral_begin() {
        return m_neutral.begin();
    }

    neutrality_iterator neutral_end() {
        return m_neutral.end();
    }

    const_neutrality_iterator neutral_begin() const {
        return m_neutral.begin();
    }

    const_neutrality_iterator neutral_end() const {
        return m_neutral.end();
    }

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

    virtual ~neutral_allele_vectorized() {}

protected:

    virtual void resize( size_t s ) {
        base_type::resize( s );

        m_neutral.resize( base_type::size() );
    }

    neutrality_vector_type  m_neutral;
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
