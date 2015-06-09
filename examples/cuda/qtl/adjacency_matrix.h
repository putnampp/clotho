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
 * DATA PAGE LAYOUT:
 *  ________________________________________________
 * | DATA (BLOCK_COUNT)                             |
 * |________________________________________________|
 * 
 * METADATA PAGE LAYOUT:
 *  ________________________________________________
 * | Free Column (1-bit per column)                 |
 * |________________________________________________|
 *
 *  ________________________________________________
 * | Free Row (1-bit per row)                       |
 * |________________________________________________|
 *
 ***********************************************************/
class adjacency_matrix : public page_manager {
public:
    static const unsigned int BIT_PER_NODE = 1;
    static const unsigned int BLOCK_PER_ROW = 32;    // (8*BLOCK_PER_ROW) -bits per row

    static const unsigned int NODE_PER_BLOCK = sizeof(page_manager::block_type) * 8 / BIT_PER_NODE; // (byte/block)*(bit/byte)/(bit/node) = (node/block)

    static const unsigned int MAX_COLUMN_NODES = NODE_PER_BLOCK * BLOCK_PER_ROW;
    static const unsigned int MAX_ROW_NODES = (page_manager::BLOCK_COUNT / BLOCK_PER_ROW);
    static const unsigned int MAX_NODE_PER_PAGE = page_manager::BLOCK_COUNT * NODE_PER_BLOCK;

    adjacency_matrix( unsigned int rows, unsigned int cols );

    size_t  free_pages_count()  const;
    size_t  row_padding()    const;
    size_t  column_padding() const;

    size_t  page_size() const;

    unsigned int column_metadata_nodes() const;
    unsigned int row_metadata_nodes() const;

    void resize( unsigned int rows, unsigned int cols );

    friend std::ostream & operator<<( std::ostream & out, const adjacency_matrix & rhs );

    virtual ~adjacency_matrix();
protected:
    unsigned int m_rows, m_cols;
    unsigned int m_row_pages, m_col_pages;
    unsigned int m_col_meta_size, m_row_meta_size;
    unsigned int m_size;

    // host
    page_type * m_metadata_pages;
    page_type * m_data_pages;

    // device
    page_type * m_device_metadata;
    page_type * m_device_data;
};

#endif  // ADJACENCY_MATRIX_H_
