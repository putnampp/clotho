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
#include "adjacency_matrix.h"




// c / (COLUMNS_PER_PAGE) = column pages
// r / (ROWS_PER_PAGE) = row pages
adjacency_matrix::adjacency_matrix( unsigned int rows, unsigned int cols ) :
    page_manager( (rows / MAX_ROW_NODES + 1) * (cols / MAX_COLUMN_NODES + 1))
    , m_rows( rows )
    , m_cols( cols )
    , m_row_pages( rows / MAX_ROW_NODES + 1)
    , m_col_pages( cols / MAX_COLUMN_NODES + 1)
{
}

size_t adjacency_matrix::free_pages_count() const {
    return 0;
}

size_t adjacency_matrix::row_padding() const {
    return m_row_pages * MAX_ROW_NODES - m_rows;
}

size_t adjacency_matrix::column_padding() const {
    return m_col_pages * MAX_COLUMN_NODES - m_cols;
}

void adjacency_matrix::resize( unsigned int rows, unsigned int cols ) {
    if( rows < m_row_pages && cols < m_col_pages ) return;

    unsigned int rpg = rows / MAX_ROW_NODES + 1;
    unsigned int cpg = cols / MAX_COLUMN_NODES + 1;

    if( rpg * cpg > capacity() ) {
        reserve( rpg * cpg );
    }

    //TODO: need to transform current matrix space

    m_row_pages = rpg;
    m_col_pages = cpg;
}

adjacency_matrix::~adjacency_matrix() {}

std::ostream & operator<<( std::ostream & out, const adjacency_matrix & rhs ) {
    out << "Total Pages: " << rhs.capacity() << " (" << rhs.allocated_size() << " bytes)" << "\n";
    out << "Page Dimensions: <" << rhs.m_row_pages << ", " << rhs.m_col_pages << ">\n";
    out << "Maximum Node Dimensions: <" << rhs.m_row_pages * adjacency_matrix::MAX_ROW_NODES << ", " << rhs.m_col_pages * adjacency_matrix::MAX_COLUMN_NODES << ">\n";
    out << "Minimum Node Dimensions: <" << rhs.m_rows << ", " << rhs.m_cols << ">\n";
    return out;
}
