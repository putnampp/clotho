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
#ifndef ADJACENCY_MATRIX_H_
#define ADJACENCY_MATRIX_H_

#include "page_manager.h"

#include <ostream>

/***********************************************************
 *
 * PAGE LAYOUT:
 * ________________________________________________
 *| Free Column bit vector ( COLUMN_HEADER_BLOCKS )|
 * ------------------------------------------------
 *| Free Rows bit vector   ( ROW_HEADER_BLOCKS )   |
 *|                                                |
 * ------------------------------------------------
 *| DATA (BLOCK_COUNT - x - y )                    |
 * ------------------------------------------------ 
 *| PADDING                                        |
 *|________________________________________________|
 *
 * BLOCK_COUNT >= COLUMN_HEADER_BLOCKS (x) + ROW_HEADER_BLOCKS (y) + NODE_PER_BLOCK * x * y
 *
 ***********************************************************/
class adjacency_matrix : public page_manager {
public:
    static const unsigned int BIT_PER_NODE = 1;
    static const unsigned int BLOCK_PER_ROW = 4;

    static const unsigned int NODE_PER_BLOCK = sizeof(page_manager::block_type) * 8 / BIT_PER_NODE; // (byte/block)*(bit/byte)/(bit/node) = (node/block)

    static const unsigned int COLUMN_HEADER_BLOCKS = BLOCK_PER_ROW;
    static const unsigned int MAX_COLUMN_NODES = NODE_PER_BLOCK * BLOCK_PER_ROW;

    static const unsigned int ROW_HEADER_BLOCKS = (page_manager::BLOCK_COUNT - COLUMN_HEADER_BLOCKS) / (NODE_PER_BLOCK * COLUMN_HEADER_BLOCKS + 1);

    static const unsigned int MAX_ROW_NODES = NODE_PER_BLOCK * ROW_HEADER_BLOCKS;
    static const unsigned int ROW_PADDING = (page_manager::BLOCK_COUNT - MAX_ROW_NODES);

    adjacency_matrix( unsigned int rows, unsigned int cols );

    size_t  free_pages_count()  const;
    size_t  row_padding()    const;
    size_t  column_padding() const;

    void resize( unsigned int rows, unsigned int cols );

    friend std::ostream & operator<<( std::ostream & out, const adjacency_matrix & rhs );

    virtual ~adjacency_matrix();
protected:
    unsigned int m_rows, m_cols;
    unsigned int m_row_pages, m_col_pages;
};

#endif  // ADJACENCY_MATRIX_H_
