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
#include "page_manager.h"

#include <cuda.h>
#include <iostream>

page_manager::page_manager( unsigned int _count) :
    m_pages(NULL)
    , m_capacity(0)
    , m_status( true )
{
    reserve( _count );
}

size_t page_manager::allocated_size() const {
    return m_capacity * PAGE_SIZE;
}

size_t page_manager::capacity() const {
    return m_capacity;
}

bool page_manager::good() const {
    return m_status;
}

void page_manager::reserve( unsigned int n ) {
    if( n == 0 || n < m_capacity ) return;

    size_t list_size = sizeof(page_type) * n; // in bytes
    if( m_pages == NULL ) {
        m_pages = (page_type *)malloc( list_size );
    } else {
        page_type * tmp_pages = m_pages;
        m_pages = (page_type *)malloc( list_size );
        memcpy( m_pages, tmp_pages, sizeof(page_type) * m_capacity);
        free(tmp_pages);
    }

    n -= m_capacity;
    while( n-- ) {
        cudaError_t err = cudaMalloc( (void **) &m_pages[m_capacity++], PAGE_SIZE );
        if( err != cudaSuccess ) {
            std::cerr << "Failed to allocate page " << (m_capacity - 1) << std::endl;
            m_status = false;
            return;
        }
    }
}

page_manager::~page_manager() {
    while( m_capacity-- ) {
        if( m_pages[m_capacity] ) {
            cudaError_t err = cudaFree( m_pages[m_capacity] );
            if( err != cudaSuccess ) {
                std::cerr << "Unsuccessfully to freed page: " << m_capacity << " (" << m_pages[m_capacity] << ")" << std::endl;
            }
        }
    }
}
