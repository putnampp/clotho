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
#ifndef SQUARE_CUDA_H_
#define SQUARE_CUDA_H_

#include <ostream>

class Square {
public:
    typedef unsigned long int_type;

    Square();

    size_t size() const;

    void operator()();

    void random_list();

    friend std::ostream & operator<<( std::ostream &, const Square & rhs );

    virtual ~Square();
protected:

    void init();
    //int_type m_a[N];

    int_type *  m_a, * m_dest;
    size_t      m_size;
    int m_maxBlocks, m_maxThreadsPerBlock;
};

#endif  // SQUARE_CUDA_H_
