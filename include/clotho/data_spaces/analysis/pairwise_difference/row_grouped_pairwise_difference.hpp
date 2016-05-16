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
#ifndef CLOTHO_ROW_GROUPED_PAIRWISE_DIFFERENCE_HPP_
#define CLOTHO_ROW_GROUPED_PAIRWISE_DIFFERENCE_HPP_

#include "clotho/data_spaces/analysis/pairwise_difference/pairwise_difference_evaluator_def.hpp"

#include "clotho/data_spaces/association_matrix/types/row_grouped.hpp"

#include "clotho/utility/popcount.hpp"

namespace clotho {
namespace genetics {

template < class BlockType, unsigned char P >
struct pairwise_difference_evaluator< BlockType, row_grouped< P > > {
    typedef BlockType   block_type;
    typedef BlockType * pointer;

    pairwise_difference_evaluator( size_t s ) {}

    size_t operator()( pointer a_start, pointer a_end, pointer b_start, pointer b_end ) {
        size_t diff = 0;
        while( a_start != a_end && b_start != b_end ) {
            block_type a = *a_start;
            block_type b = *b_start;

            diff += popcount( a ^ b );

            a_start += row_grouped< P >::GROUP_SIZE;
            b_start += row_grouped< P >::GROUP_SIZE;
        }

        assert( a_start == a_end && b_start == b_end );

        return diff;
    }
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_ROW_GROUPED_PAIRWISE_DIFFERENCE_HPP_
