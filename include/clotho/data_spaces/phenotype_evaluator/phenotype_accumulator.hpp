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

template < class TraitSpaceType >
class phenotype_accumulator {
public:
    typedef phenotype_accumulator< TraitSpaceType >     self_type;
    typedef TraitSpaceType                              trait_space_type;
    typedef typename trait_space_type::weight_type      weight_type;

    typedef typename trait_space_type::const_iterator   const_iterator;
/**
 *
 * rows => alleles
 * columns => traits
 */
    phenotype_accumulator( trait_space_type * traits,  weight_type * results ) :
        m_traits( traits )
        , m_results( results )
    {}

    phenotype_accumulator( const self_type & other ) :
        m_traits( other.m_traits )
        , m_results( other.m_results )
    {}

    void operator()( unsigned int row ) {
        const_iterator it = m_traits->begin( row );
        const_iterator end = m_traits->end( row );

        weight_type * tmp = m_results;
        while( it != end ) {
            (*tmp++) += (*it++);
        }
    }

    weight_type * getResults() {
        return m_results;
    }

    virtual ~phenotype_accumulator() {}

protected:
    trait_space_type    * m_traits;
    weight_type         * m_results;
};

}   // namespace genetics
}   // namespace clotho
#endif  // CLOTHO_PHENOTYPE_ACCUMULATOR_HPP_

