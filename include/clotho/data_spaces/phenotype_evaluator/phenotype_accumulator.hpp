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
#ifndef CLOTHO_PHENOTYPE_ACCUMULATOR_HPP_
#define CLOTHO_PHENOTYPE_ACCUMULATOR_HPP_

namespace clotho {
namespace genetics {

template < class WeightType >
class phenotype_accumulator {
public:
    typedef phenotype_accumulator< WeightType > self_type;
    typedef WeightType  weight_type;

/**
 *
 * rows => alleles
 * columns => traits
 */
    phenotype_accumulator( weight_type * weights, unsigned int rows, unsigned int columns, weight_type * results ) :
        m_weights( weights )
        , m_rows(rows)
        , m_columns( columns )
        , m_results( results )
    {}

    phenotype_accumulator( const self_type & other ) :
        m_weights( other.m_weights )
        , m_rows( other.m_rows )
        , m_columns( other.m_columns )
        , m_results( other.m_results )
    {}

    void operator()( unsigned int row ) {
#ifdef DEBUGGING
        assert( row < m_rows );
#endif  // DEBUGGING

        row *= m_columns;
        for( unsigned int i = 0; i < m_columns; ++i ) {
            m_results[ i ] += m_weights[ row + i ];
        }
    }

    virtual ~phenotype_accumulator() {}

protected:
    weight_type     * m_weights;

    unsigned int    m_rows, m_columns;

    weight_type     * m_results;
};

}   // namespace genetics
}   // namespace clotho
#endif  // CLOTHO_PHENOTYPE_ACCUMULATOR_HPP_

