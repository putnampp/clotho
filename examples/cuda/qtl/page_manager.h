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
#ifndef PAGE_MANAGER_H_
#define PAGE_MANAGER_H_

#include <cuda_runtime.h>

class page_manager {
public:
    typedef unsigned int block_type;
    typedef block_type * page_type;

    static const unsigned int PAGE_SIZE = 4096;
    static const unsigned int BLOCK_COUNT = PAGE_SIZE / sizeof(block_type);
    

    page_manager( unsigned int _count = 1 );

/**
 * Returns the total number of bytes allocated
 */
    size_t  allocated_size()    const;

/**
 * Returns the number of allocated pages
 */
    size_t  capacity()          const;

    bool good() const;

    void reserve( unsigned int n );

    virtual ~page_manager();
protected:
    page_type * m_pages;
    unsigned int m_capacity;
    bool m_status;
};
#endif  // PAGE_MANAGER_H_
