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
#ifndef CLOTHO_COLUMN_ALIGNMENT_PAIRWISE_DIFFERENCE_HPP_
#define CLOTHO_COLUMN_ALIGNMENT_PAIRWISE_DIFFERENCE_HPP_

#include "clotho/data_spaces/analysis/pairwise_difference/pairwise_difference_evaluator_def.hpp"

#include "clotho/data_spaces/association_matrix/types/column_aligned.hpp"

#include "clotho/utility/popcount.hpp"

namespace clotho {
namespace genetics {

template < class BlockType >
struct pairwise_difference_evaluator< BlockType, column_aligned > {
    typedef BlockType block_type;
    typedef BlockType * pointer;

    size_t m_step;

    pairwise_difference_evaluator( size_t s ) : m_step(s) {}

    size_t operator()( pointer a_start, pointer a_end, pointer b_start, pointer b_end ) {
        size_t diff = 0;
        while( a_start != a_end && b_start != b_end ) {
            block_type a = *a_start;
            block_type b = *b_start;

            diff += popcount( a ^ b );

            a_start += m_step;
            b_start += m_step;
        }

        assert( a_start == a_end && b_start == b_end );
        return diff;
    }
};

}   // namespace genetics
}   // namespace clotho
#endif  // CLOTHO_COLUMN_ALIGNMENT_PAIRWISE_DIFFERENCE_HPP_
