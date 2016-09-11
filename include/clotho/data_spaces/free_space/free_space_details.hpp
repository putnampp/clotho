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
#ifndef CLOTHO_FREE_SPACE_DETAILS_HPP_
#define CLOTHO_FREE_SPACE_DETAILS_HPP_

#include "clotho/utility/bit_helper.hpp"
#include "clotho/utility/debruijn_bit_walker.hpp"

namespace clotho {
namespace genetics {

template < class SizeType >
class free_space_details {
public:
    typedef SizeType                                        size_type;

    typedef size_type *                                     index_vector;
    typedef size_type *                                     iterator;
    typedef size_type *                                     const_iterator;

    free_space_details( ) :
        m_indices( NULL )
        , m_fixed_count(0)
        , m_lost_count(0)
        , m_free_count(0)
        , m_width(0)
        , m_size(0)
    {}

    size_type variable_count() const {
        return m_width - m_free_count;
    }

    size_type free_size() const {
        return m_free_count;
    }

    iterator free_begin() {
        return m_indices;
    }

    iterator free_end() {
        return m_indices + m_free_count;
    }

    const_iterator free_begin() const {
        return m_indices;
    }

    const_iterator free_end() const {
        return m_indices + m_free_count;
    }

    size_type fixed_size() const {
        return m_fixed_count;
    }

    iterator fixed_begin() {
        return m_indices + m_width;
    }

    iterator fixed_end() {
        return m_indices + m_width + m_fixed_count;
    }

    const_iterator fixed_begin() const {
        return m_indices + m_width;
    }

    const_iterator fixed_end() const {
        return m_indices + m_width + m_fixed_count;
    }

    size_type lost_size() const {
        return m_lost_count;
    }

    iterator lost_begin() {
        return m_indices + 2 * m_width;
    }

    iterator lost_end() {
        return m_indices + 2 * m_width + m_lost_count;
    }

    const_iterator lost_begin() const {
        return m_indices + 2 * m_width;
    }

    const_iterator lost_end() const {
        return m_indices + 2 * m_width + m_lost_count;
    }

    virtual ~free_space_details() {
        if( m_indices != NULL ) {
            delete [] m_indices;
        }
    }

protected:
    
    void resize( size_type s ) {
        if( s > m_width ) {
            if( m_indices != NULL ) {
                delete [] m_indices;
            }

            m_indices = new size_type[ 3 * s ];
            m_size = 3 * s;
            m_width = s;
        }
    }
    template < class BlockType >
    void analyze_free_indices( BlockType * fixed, BlockType * variable, unsigned int W, unsigned int M ) {
        typedef BlockType                                           block_type;
        typedef clotho::utility::debruijn_bit_walker< block_type >  bit_walker_type;
        typedef clotho::utility::BitHelper< block_type >            bit_helper_type;

        size_type * _fix = m_indices + m_width;
        size_type * _lost = _fix + m_width; 

        m_free_count = 0;
        m_lost_count = 0;
        m_fixed_count = 0;

        unsigned int j = 0;
        for( unsigned int i = 0; i < W; ++i ) {
            block_type fx = fixed[ i ];
            block_type var = variable[ i ];

            block_type ls = ~(fx | var);

            while( fx ) {
                size_type b_idx = bit_walker_type::unset_next_index( fx ) + j;
                if( b_idx < M ) {
                    m_indices[ m_free_count++ ] = b_idx;
                    _fix[ m_fixed_count++ ] = b_idx;
                }
            }

            while( ls ) {
                size_type idx = bit_walker_type::unset_next_index( ls ) + j;
                if( idx < M ) {
                    m_indices[ m_free_count++ ] = idx;
                    _lost[ m_lost_count++ ] = idx;
                }
            }

            j += bit_helper_type::BITS_PER_BLOCK;
        }
    }

    size_type  * m_indices;
    size_type  m_fixed_count, m_lost_count, m_free_count;
    size_type  m_width, m_size;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_FREE_SPACE_DETAILS_HPP_
