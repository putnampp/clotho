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
#ifndef CLOTHO_ROW_GROUPED_BLOCK_INDEXING_DEF_HPP_
#define CLOTHO_ROW_GROUPED_BLOCK_INDEXING_DEF_HPP_

#include "clotho/data_spaces/association_matrix/indexing/block_indexing_def.hpp"

#include "clotho/data_spaces/association_matrix/types/row_grouped.hpp"

namespace clotho {
namespace genetics {

template < class BlockType, unsigned char P >
struct block_indexing< BlockType, row_grouped< P > > {

/**
 * bpr - blocks_per_row
 */
    inline size_t block_index( size_t row, size_t col, size_t bpr ) const {
        return block_row_offset( row, bpr ) + block_column_offset( col );
    }

    inline size_t block_row_offset( size_t row, size_t bpr ) const {
        return (row / row_grouped< P >::GROUP_SIZE) * bpr * row_grouped< P >::GROUP_SIZE + (row % row_grouped< P >::GROUP_SIZE);
    }

    inline size_t block_column_offset( size_t col ) const {
        return (col / clotho::utility::BitHelper< BlockType >::BITS_PER_BLOCK) * row_grouped< P >::GROUP_SIZE;
    }
};

}   // namespace genetics
}   // namespace clotho
#endif  // CLOTHO_ROW_GROUPED_BLOCK_INDEXING_DEF_HPP_
