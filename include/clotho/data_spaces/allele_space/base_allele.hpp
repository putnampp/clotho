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
#ifndef CLOTHO_BASE_ALLELE_HPP_
#define CLOTHO_BASE_ALLELE_HPP_

#include <vector>

#include "clotho/data_spaces/allele_space/growable1D.hpp"

#include <iostream>

namespace clotho {
namespace genetics {

template < class PositionType = double >
class base_allele_vectorized : public growable1D {
public:
    typedef PositionType     position_type;

    typedef std::vector< position_type >                     position_vector_type;
    typedef typename position_vector_type::iterator          position_iterator;
    typedef typename position_vector_type::const_iterator    const_position_iterator;

    position_type getPositionAt( size_t index ) const {
        return m_positions[ index ];
    }

    void setPositionAt( size_t index, position_type pos ) {
        if( index >= m_positions.size() ) {
            resize( index + 1 );
        }

        m_positions[ index ] = pos;
    }

    position_iterator position_begin() {
        return m_positions.begin();
    }

    position_iterator position_end() {
        return m_positions.end();
    }

    const_position_iterator position_begin() const {
        return m_positions.begin();
    }

    const_position_iterator position_end() const {
        return m_positions.begin();
    }

    virtual size_t grow( size_t s ) {
        resize( s );

        return size();
    }

    size_t size() const {
        return m_positions.size();
    }
    
    virtual ~base_allele_vectorized() {}

protected:

    virtual void resize( size_t s ) {
        std::cout << "base 1D resize" << std::endl;
        m_positions.resize( s );
    }

    position_vector_type         m_positions;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_BASE_ALLELE_HPP_
