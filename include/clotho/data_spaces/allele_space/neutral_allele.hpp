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

namespace clotho {
namespace genetics {

template < class PositionType >
class neutral_allele_vectorized : public base_allele_vectorized< PositionType > {
public:

    typedef base_allele_vectorized< PositionType >   base_type;

    typedef std::vector< bool >     neutrality_vector_type;
    typedef typename neutrality_vector_type::iterator       neutrality_iterator;
    typedef typename neutrality_vector_type::const_iterator const_neutrality_iterator;

    bool getNeutralAt( size_t index ) const {
        return m_neutral[ index ];
    }

    void setNeutralAt( size_t index, bool neu ) {
        if( index >= base_type::size() ) {
            resize( index + 1 );
        }

        m_neutral[ index ] = neu;
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

#endif  // CLOTHO_NEUTRAL_ALLELE_HPP_
