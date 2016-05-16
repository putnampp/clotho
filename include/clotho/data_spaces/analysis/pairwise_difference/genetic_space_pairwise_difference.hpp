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
#ifndef GENETIC_SPACE_PAIRWISE_DIFFERENCE_HPP_
#define GENETIC_SPACE_PAIRWISE_DIFFERENCE_HPP_

#include "clotho/data_spaces/analysis/pairwise_difference/association_matrix_pairwise_difference.hpp"
#include "clotho/data_spaces/population_space/genetic_space.hpp"

namespace clotho {
namespace genetics {

template < class AlleleType, class BlockType, class OrganizationType >
class pairwise_difference< genetic_space< AlleleType, BlockType, OrganizationType > > : public pairwise_difference< typename genetic_space< AlleleType, BlockType, OrganizationType >::association_type > {
public:
    typedef pairwise_difference< typename genetic_space< AlleleType, BlockType, OrganizationType >::association_type >    base_type;
    typedef genetic_space< AlleleType, BlockType, OrganizationType >                                       space_type;

    void evaluate( space_type & ss ) {
        base_type::evaluate( ss.getSequenceSpace() );
    }

    template < class Iterator >
    void evaluate( space_type & ss, Iterator first, Iterator last ) {
        base_type::evaluate( ss.getSequenceSpace(), first, last );
    }

    virtual ~pairwise_difference() {}
};

}   // namespace genetics
}   // namespace clotho

#endif  // GENETIC_SPACE_PAIRWISE_DIFFERENCE_HPP_
