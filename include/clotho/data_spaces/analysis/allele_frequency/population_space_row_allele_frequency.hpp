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
#ifndef POPULATION_SPACE_ROW_ALLELE_FREQUENCY_HPP_
#define POPULATION_SPACE_ROW_ALLELE_FREQUENCY_HPP_

#include <boost/property_tree/ptree.hpp>
#include "clotho/data_spaces/population_space/population_space_row.hpp"

#include "clotho/data_spaces/analysis/allele_frequency/margin_details.hpp"
#include "clotho/data_spaces/analysis/allele_frequency/frequency_evaluator.hpp"

namespace clotho {
namespace genetics {

template < class BlockType, class WeightType, class SizeType >
class allele_frequency< population_space_row< BlockType, WeightType >, SizeType > : public margin_details< SizeType >  {
public:
    typedef population_space_row< BlockType, WeightType > space_type;
    typedef typename margin_details< SizeType >::size_type size_type;
    typedef frequency_evaluator< space_type, size_type > evaluator_type;

    allele_frequency() {}

    void evaluate( space_type & ss ) {
        this->resize( ss.getMaxAlleles(), ss.haploid_genome_count() );

        this->m_indices.reset();
        this->m_indices.flip();

        evaluator_type eval;

        eval( ss, this->m_indices, this->m_column_margin, this->m_row_margin );
    }

/**
 * Iterator sub population
 */
    template < class Iterator >
    void evaluate( space_type & ss, Iterator first, Iterator last ) {
        this->resize( ss.getMaxAlleles(), ss.haploid_genome_count() );

        this->m_indices.reset();

        size_type N = 0;
        while( first != last ) {
            size_type i = *first++;
            this->m_indices[i] = true;
            ++N;
        }

        evaluator_type eval;

        eval( ss, this->m_indices, this->m_column_margin, this->m_row_margin );
    }

    virtual ~allele_frequency() {} 
};

}   // namespace genetics
}   // namespace clotho
#endif  // POPULATION_SPACE_ROW_ALLELE_FREQUENCY_HPP_
