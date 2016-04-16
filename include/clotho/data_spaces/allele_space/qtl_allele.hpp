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
#ifndef CLOTHO_QTL_ALLELE_HPP_
#define CLOTHO_QTL_ALLELE_HPP_

#include "clotho/data_spaces/allele_space/neutral_allele.hpp"
#include "clotho/data_spaces/growable2D.hpp"

#include <iostream>

namespace clotho {
namespace genetics {

template < class PositionType, class WeightType >
class qtl_allele_vectorized : public neutral_allele_vectorized< PositionType >, public growable2D {
public:
    typedef neutral_allele_vectorized< PositionType > base_type;

    typedef WeightType  weight_type;

    qtl_allele_vectorized( size_t rows = 1, size_t columns = 1 ) :
        m_weights(NULL)
        , m_rows(0)
        , m_columns(0)
        , m_size(0)
    {
        this->resize(rows, columns);
    }
    

    virtual size_t grow( size_t alleles ) {
        this->resize( alleles );
        return m_rows;
    }

    virtual size_t grow( size_t alleles, size_t traits ) {
        this->resize( alleles, traits );

        return m_rows;
    }

    size_t allele_count() const {
        return m_rows;
    }

    size_t trait_count() const {
        return m_columns;
    }

    size_t allocated_size() const {
        return m_size;
    }

    virtual ~qtl_allele_vectorized() {
        if( m_weights != NULL ) {
            delete [] m_weights;
        }
    }

protected:

    virtual void resize( size_t rows, size_t columns ) {
        std::cout << "QTL 2D resize" << std::endl;
        base_type::resize( rows );

        rows = base_type::size();

        size_t new_size = rows * columns;

        assert( 0 < new_size );

        if( m_size < new_size ) {
            if( m_weights != NULL ) {
                delete [] m_weights;
            }

            m_weights = new weight_type[ new_size ];

            m_size = new_size;
        }

        m_columns = columns;
        m_rows = rows;
    }

    virtual void resize( size_t rows ) {
        std::cout << "QTL 1D resize" << std::endl;    
        this->resize( rows, m_columns );
    }

    weight_type * m_weights;
    size_t      m_rows, m_columns, m_size;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_QTL_ALLELE_HPP_
