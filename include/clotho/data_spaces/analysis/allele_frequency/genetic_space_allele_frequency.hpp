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
#ifndef GENETIC_SPACE_ALLELE_FREQUENCY_HPP_
#define GENETIC_SPACE_ALLELE_FREQUENCY_HPP_

#include "clotho/data_spaces/analysis/allele_frequency/association_matrix_allele_frequency.hpp"
#include "clotho/data_spaces/population_space/genetic_space.hpp"

namespace clotho {
namespace genetics {

template < class AlleleType, class BlockType, class AlignmentType >
class allele_frequency< genetic_space< AlleleType, BlockType, AlignmentType > > : public allele_frequency< typename genetic_space< AlleleType, BlockType, AlignmentType >::association_type > {
public:
    typedef genetic_space< AlleleType, BlockType, AlignmentType >  genetic_space_type;
    typedef allele_frequency< typename genetic_space< AlleleType, BlockType, AlignmentType >::association_type > base_type;

    void evaluate( genetic_space_type & gs ) {
        base_type::evaluate( gs.getSequenceSpace() );
    }

    template < class Iterator >
    void evaluate( genetic_space_type & gs, Iterator first, Iterator last ) {
        base_type::evaluate( gs.getSequenceSpace(), first, last );
    }

    virtual ~allele_frequency() {}
};

}   // namespace genetics
}   // namespace clotho

#endif  // GENETIC_SPACE_ALLELE_FREQUENCY_HPP_
