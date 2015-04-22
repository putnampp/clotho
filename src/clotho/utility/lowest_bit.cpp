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
#include "clotho/utility/lowest_bit.h"

namespace clotho {
namespace utility {

lowest_bit_256::lowest_bit_256() {
    initialize();
}

void lowest_bit_256::initialize() {
    block_type val = (block_type)( max_values - 1);
    while( val ) {
        block_type idx = 0, tmp = val;
        while( tmp ) {
            if( tmp & 1 ) {
                m_vals[ val ].bit_index = idx++;
                tmp >>= 1;
                break;
            }

            tmp >>= 1;
            ++idx;
        }

        if( tmp ) {
            do {
                if( tmp & 1 ) {
                    m_vals[ val ].bit_shift_next = idx;
                    m_vals[ val ].next_ptr = &m_vals[tmp];
                    break;
                }
                tmp >>= 1;
                ++idx;
            } while( tmp );
        } else {
            m_vals[ val ].bit_shift_next = idx;
            m_vals[ val ].next_ptr = NULL;
        }
        --val;
    }
}

const lowest_bit_256::value_type & lowest_bit_256::operator[]( block_type val ) const {
    return m_vals[ val ];
}

unsigned char lowest_bit_256::low_bit_index( block_type val ) const {
    return m_vals[ val ].bit_index;
}

unsigned char lowest_bit_256::next_shift_offset( block_type val ) const {
    return m_vals[ val ].bit_shift_next;
}

const lowest_bit_256::value_type * lowest_bit_256::begin() const {
    return m_vals;
}

lowest_bit_256::~lowest_bit_256() {}

lowest_bit_64K::lowest_bit_64K() {
    initialize();
}

void lowest_bit_64K::initialize() {
    block_type val = (max_values - 1);
    while( val ) {
        block_type idx = 0, tmp = val;
        while( tmp ) {
            if( tmp & 1 ) {
                m_vals[ val ].bit_index = idx++;
                tmp >>= 1;
                break;
            }

            tmp >>= 1;
            ++idx;
        }

        if( tmp ) {
            do {
                if( tmp & 1 ) {
                    m_vals[ val ].bit_shift_next = idx;
                    m_vals[ val ].next_ptr = &m_vals[tmp];
                    break;
                }
                tmp >>= 1;
                ++idx;
            } while( tmp );
        } else {
            m_vals[ val ].bit_shift_next = idx;
            m_vals[ val ].next_ptr = NULL;
        }
        --val;
    }
}

const lowest_bit_64K::value_type & lowest_bit_64K::operator[]( block_type val ) const {
    return m_vals[ val ];
}

unsigned char lowest_bit_64K::low_bit_index( block_type val ) const {
    return m_vals[ val ].bit_index;
}

unsigned char lowest_bit_64K::next_shift_offset( block_type val ) const {
    return m_vals[ val ].bit_shift_next;
}

const lowest_bit_64K::value_type * lowest_bit_64K::begin() const {
    return m_vals;
}

lowest_bit_64K::~lowest_bit_64K() {}

}   // namespace utility
}   // namespace clotho
