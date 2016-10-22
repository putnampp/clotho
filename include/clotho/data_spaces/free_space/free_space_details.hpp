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

/*
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

        size_type * _free = m_indices;
        size_type * _fix = m_indices + m_width;
        size_type * _lost = _fix + m_width; 

        for( unsigned int i = 0; i < W; ++i ) {
            block_type fx = fixed[ i ];
            block_type var = variable[ i ];

            block_type ls = ~(fx | var);

            size_type idx = i * bit_helper_type::BITS_PER_BLOCK;
            while( fx ) {
                idx += bit_walker_type::next_and_shift( fx );
                if( idx < M ) {
                    *_free++ = idx;
                    *_fix++ = idx;
                }
            }

            idx = i * bit_helper_type::BITS_PER_BLOCK;
            while( ls ) {
                idx += bit_walker_type::next_and_shift( ls );
                if( idx < M ) {
                    *_free++ = idx;
                    *_lost++ = idx;
                }
            }
        }

        m_free_count = (_free - m_indices);
        m_fixed_count = (_fix - (m_indices + m_width));
        m_lost_count = (_lost - (m_indices + 2 * m_width));

#ifdef DEBUGGING
        BOOST_LOG_TRIVIAL(debug) << "Free count: " << m_free_count << "; Fixed count: " << m_fixed_count << "; Lost count: " << m_lost_count;
#endif  // DEBUGGING
    }

    size_type  * m_indices;
    size_type  m_fixed_count, m_lost_count, m_free_count;
    size_type  m_width, m_size;
};
*/

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
        , m_size(0)
        , m_allocated_size(0)
    {}

    size_type variable_count() const {
        return m_size - (m_fixed_count + m_lost_count);
    }

    size_type free_size() const {
        return m_fixed_count + m_lost_count;
    }

    iterator free_begin() {
        return m_indices;
    }

    iterator free_end() {
        return m_indices + m_fixed_count + m_lost_count;
    }

    const_iterator free_begin() const {
        return m_indices;
    }

    const_iterator free_end() const {
        return m_indices + m_fixed_count + m_lost_count;
    }

    size_type fixed_size() const {
        return m_fixed_count;
    }

    iterator fixed_begin() {
        return m_indices;
    }

    iterator fixed_end() {
        return m_indices + m_fixed_count;
    }

    const_iterator fixed_begin() const {
        return m_indices;
    }

    const_iterator fixed_end() const {
        return m_indices + m_fixed_count;
    }

    size_type lost_size() const {
        return m_lost_count;
    }

    iterator lost_begin() {
        return m_indices + m_fixed_count;
    }

    iterator lost_end() {
        return m_indices + m_fixed_count + m_lost_count;
    }

    const_iterator lost_begin() const {
        return m_indices + m_fixed_count;
    }

    const_iterator lost_end() const {
        return m_indices + m_fixed_count + m_lost_count;
    }

    virtual ~free_space_details() {
        if( m_indices != NULL ) {
            delete [] m_indices;
        }
    }

protected:
    
    void resize( size_type s ) {
        if( s > m_allocated_size ) {
            if( m_indices != NULL ) {
                delete [] m_indices;
            }

            m_indices = new size_type[ s ];
            m_allocated_size = s;
        }
        m_size = s;
    }

    template < class BlockType >
    void analyze_free_indices( BlockType * fixed_var, BlockType * end, unsigned int M ) {
        typedef BlockType                                           block_type;
        typedef clotho::utility::debruijn_bit_walker< block_type >  bit_walker_type;
        typedef clotho::utility::BitHelper< block_type >            bit_helper_type;

        block_type * tmp = fixed_var;
        unsigned int idx = 0;
        size_type offset = 0;
        while( tmp != end ) {
            block_type fx = *tmp;

            unsigned int j = idx;
            while( fx ) {
                j += bit_walker_type::next_and_shift( fx );
#ifdef DEBUGGING
                // do not push sequence padding bits
                assert( j < M );
#endif  // DEBUGGING

                m_indices[ offset++ ] = j;
            }

            idx += bit_helper_type::BITS_PER_BLOCK;
            tmp += 2;
        }

        m_fixed_count = offset;

        tmp = fixed_var + 1;
        end += 1;
        idx = 0;
        while( tmp != end ) {
            block_type ls = *tmp;

            ls = ~ls;
            unsigned int j = idx;
            while( ls ) {
                j += bit_walker_type::next_and_shift( ls );
#ifdef DEBUGGING
                // do not push sequence padding bits
                assert( j < M );
#endif  // DEBUGGING

                m_indices[ offset++ ] = j;
            }

            idx += bit_helper_type::BITS_PER_BLOCK;
            tmp += 2;
        }

        m_lost_count = offset - m_fixed_count;

#ifdef DEBUGGING
        std::cerr << "lost_count = " << offset << " - " << m_fixed_count << " = " << m_lost_count << std::endl;
        BOOST_LOG_TRIVIAL(debug) << "Free count: " << free_size() << "; Fixed count: " << fixed_size() << "; Lost count: " << lost_size();
        std::cerr << "Free count: " << free_size() << "; Fixed count: " << fixed_size() << "; Lost count: " << lost_size() << std::endl;
#endif  // DEBUGGING
    }

    template < class BlockType >
    void analyze_free_indices( BlockType * fixed, BlockType * variable, unsigned int W, unsigned int M ) {
        typedef BlockType                                           block_type;
        typedef clotho::utility::debruijn_bit_walker< block_type >  bit_walker_type;
        typedef clotho::utility::BitHelper< block_type >            bit_helper_type;

        block_type * first = fixed, * last = fixed + W;

        unsigned int offset = 0;

        unsigned int idx = 0;
        while( first != last ) {
            block_type fx = *first++;

            unsigned int b_idx = idx;
            while( fx ) {
                b_idx += bit_walker_type::next_and_shift( fx );
                if( b_idx < M ) {
                    m_indices[ offset++ ] = b_idx;
                }
            }
            idx += bit_helper_type::BITS_PER_BLOCK;
        }

        m_fixed_count = offset;

        first = variable;
        last = variable + W;
        idx = 0;
        while( first != last ) {
            block_type ls = *first++;
            ls = ~ls;

            unsigned int b_idx = idx;
            while( ls ) {
                b_idx += bit_walker_type::next_and_shift( ls );
                if( b_idx < M ) {
                    m_indices[ offset++ ] = b_idx;
                }
            }

            idx += bit_helper_type::BITS_PER_BLOCK;
        }

        m_lost_count = offset - m_fixed_count;

#ifdef DEBUGGING
        std::cerr << "lost_count = " << offset << " - " << m_fixed_count << " = " << m_lost_count << std::endl;
        BOOST_LOG_TRIVIAL(debug) << "Free count: " << free_size() << "; Fixed count: " << fixed_size() << "; Lost count: " << lost_size();
        std::cerr << "Free count: " << free_size() << "; Fixed count: " << fixed_size() << "; Lost count: " << lost_size() << std::endl;
#endif  // DEBUGGING
    }

    size_type  * m_indices;
    size_type  m_fixed_count, m_lost_count;
    size_type  m_size, m_allocated_size;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_FREE_SPACE_DETAILS_HPP_
