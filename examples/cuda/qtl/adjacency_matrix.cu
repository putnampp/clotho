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

#include <iostream>

adjacency_matrix::adjacency_matrix( unsigned int rows, unsigned int cols ) :
    page_manager( 0 )
    , m_rows( 0 )
    , m_cols( 0 )
    , m_row_pages( 0 )
    , m_col_pages( 0 )
    , m_col_meta_size(0)
    , m_row_meta_size(0)
    , m_size(0)
    , m_metadata_pages(NULL)
    , m_data_pages(NULL)
{
    resize( rows, cols );
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

size_t adjacency_matrix::page_size() const {
    return m_col_meta_size + m_row_meta_size + m_row_pages * m_col_pages;
}

void adjacency_matrix::resize( unsigned int rows, unsigned int cols ) {
    if( rows < m_row_pages && cols < m_col_pages ) return;

    unsigned int rpg = rows / MAX_ROW_NODES + (( rows % MAX_ROW_NODES ) ? 1 : 0 );
    unsigned int cpg = cols / MAX_COLUMN_NODES + (( cols % MAX_COLUMN_NODES ) ? 1 : 0 );

    unsigned int meta_cols = cpg * MAX_COLUMN_NODES;
    unsigned int meta_rows = rpg * MAX_ROW_NODES;

    unsigned int mcpg = meta_cols / MAX_NODE_PER_PAGE + ((meta_cols % MAX_NODE_PER_PAGE ) ? 1 : 0);
    unsigned int mrpg = meta_rows / MAX_NODE_PER_PAGE + ((meta_rows % MAX_NODE_PER_PAGE ) ? 1 : 0);

    unsigned int dpg = rpg*cpg;
    unsigned int mpg = mcpg + mrpg;
    unsigned int totpages = mpg + dpg;

    reserve( totpages );

    //TODO: need to transform current matrix space
    // Current algorithm re-maps all pages

    if( mpg > m_col_meta_size + m_row_meta_size ) {
        if( m_metadata_pages ) {
            free( m_metadata_pages );

            cudaError_t err = cudaFree(m_device_metadata);
            if( err != cudaSuccess ) {
                std::cerr << "Unable to free device metadata\n";
            }
        }
        m_metadata_pages = (page_type *) malloc( sizeof(page_type) * mpg);
        cudaError_t err = cudaMalloc(&m_device_metadata, sizeof(page_type)*mpg);
        if( err != cudaSuccess ) {
            std::cerr << "Unable to allocate device metadata\n";
        }
    }

    if( dpg > m_row_pages * m_col_pages ) {
        if( m_data_pages ) {
            free( m_data_pages );
            
            cudaError_t err =  cudaFree( m_device_data);
            if( err != cudaSuccess ) {
                std::cerr << "Unable to free device data\n";
            }
        }

        m_data_pages = (page_type *)malloc( sizeof(page_type) * dpg);
        cudaError_t err = cudaMalloc(&m_device_data, sizeof(page_type) * dpg);
        if( err != cudaSuccess ) {
            std::cerr << "Unable to allocate device data\n";
        }
    }

    m_size = 0;
    m_col_meta_size = mcpg;
    m_row_meta_size = mrpg;
    memcpy( m_metadata_pages, m_pages + m_size, sizeof(page_type) * mpg);
    cudaMemcpy(m_device_metadata, m_metadata_pages, sizeof(page_type) * mpg, cudaMemcpyHostToDevice);

    m_size += mpg;
    m_row_pages = rpg;
    m_col_pages = cpg;
    memcpy(m_data_pages, m_pages + m_size, sizeof(page_type) * dpg);
    cudaMemcpy(m_device_data, m_data_pages, sizeof(page_type) * dpg, cudaMemcpyHostToDevice);
    m_size += dpg;

    m_rows = rows;
    m_cols = cols;
}

unsigned int adjacency_matrix::column_metadata_nodes() const {
    return m_col_meta_size * MAX_NODE_PER_PAGE;
}

unsigned int adjacency_matrix::row_metadata_nodes() const {
    return m_row_meta_size * MAX_NODE_PER_PAGE;
}

adjacency_matrix::~adjacency_matrix() {
    if( m_metadata_pages ) {
        free( m_metadata_pages );
        cudaFree( m_device_metadata );
    }

    if( m_data_pages ) {
        free( m_data_pages );
        cudaFree( m_device_data );       
    }
}

std::ostream & operator<<( std::ostream & out, const adjacency_matrix & rhs ) {
    out << "Total Pages: " << rhs.capacity() << " (" << rhs.allocated_size() << " bytes)" << "\n";
    out << "Data Page Dimensions: <" << rhs.m_row_pages << ", " << rhs.m_col_pages << ">\n";
    out << "Maximum Node Dimensions: <" << rhs.m_row_pages * adjacency_matrix::MAX_ROW_NODES << ", " << rhs.m_col_pages * adjacency_matrix::MAX_COLUMN_NODES << ">\n";
    out << "Minimum Node Dimensions: <" << rhs.m_rows << ", " << rhs.m_cols << ">\n";

    out << "Column Metadata Nodes: " << rhs.column_metadata_nodes() << "\n";
    out << "Row Metadata Nodes: " << rhs.row_metadata_nodes() << "\n";
    return out;
}
