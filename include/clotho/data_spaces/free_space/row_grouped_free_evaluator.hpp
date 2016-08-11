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
#ifndef CLOTHO_ROW_GROUPED_FREE_SPACE_EVALUATOR_HPP_
#define CLOTHO_ROW_GROUPED_FREE_SPACE_EVALUATOR_HPP_

#include "clotho/data_spaces/association_matrix/row_grouped_association_matrix.hpp"
#include "clotho/utility/bit_helper.hpp"
#include "clotho/utility/debruijn_bit_walker.hpp"

#include <algorithm>

namespace clotho {
namespace genetics {

template < class BlockType, unsigned char P >
class free_space_evaluator< association_matrix< BlockType, row_grouped< P > > > {
public:
    typedef association_matrix< BlockType, row_grouped< P > > space_type;
    typedef typename space_type::block_type                 block_type;
    typedef size_t *                                        result_type;

    typedef typename space_type::raw_block_pointer    block_pointer;

    typedef clotho::utility::debruijn_bit_walker< block_type >  bit_walker_type;
    typedef clotho::utility::BitHelper< block_type >            bit_helper_type;

    void operator()( space_type & ss, result_type res, size_t & fixed_offset, size_t & lost_offset, size_t & free_count, size_t M ) {

        size_t W = bit_helper_type::padded_block_count( M );

        block_type * tmp = new block_type[ 2 * W ];

        memset( tmp, 255, sizeof(block_type) * W );
        memset( tmp + W, 0, sizeof(block_type) * W );

        const size_t row_count = ss.block_row_count();

        size_t fo = fixed_offset;
        size_t lo = lost_offset;
        size_t fr = free_count;
        
        block_type * fx_ptr = tmp;
        block_type * var_ptr = tmp + W;

        size_t k = 0;
        while( k < row_count ) {

            block_pointer start = ss.begin_block_row( k );
            block_pointer end = ss.end_block_row( k );

            fx_ptr = tmp;
            var_ptr = tmp + W;

            size_t i = 0; 
            while( start != end ) {
                block_type b = *start++;
                *fx_ptr &= b;
                *var_ptr |= b;
                if( ++i == row_grouped< P >::GROUP_SIZE ) {
                    ++fx_ptr;
                    ++var_ptr;
                    i = 0;
                }
            }

            assert(fx_ptr == tmp + W);

            k += row_grouped< P >::GROUP_SIZE;
        }

        size_t j = 0;
        
        fx_ptr = tmp;
        var_ptr = tmp + W;
        while( fx_ptr != tmp + W ) {
            block_type fx = *fx_ptr++;
            block_type var = *var_ptr++;

            block_type ls = ~(fx | var);

            while( fx ) {
                size_t b_idx = bit_walker_type::unset_next_index( fx ) + j;
                if( b_idx < M ) {
                    res[ fr++ ] = b_idx;
                    res[ fo++ ] = b_idx;
                }
            }

            while( ls ) {
                size_t idx = bit_walker_type::unset_next_index( ls ) + j;
                if( idx < M ) {
                    res[ fr++ ] = idx;
                    res[ lo++ ] = idx;
                }
            }

            j += bit_helper_type::BITS_PER_BLOCK;
        }

        fixed_offset = fo;
        lost_offset = lo;
        free_count = fr;

        delete [] tmp;
    }
};

template < class BlockType >
class free_space_evaluator< association_matrix< BlockType, row_grouped< 1 > > > {
public:
    typedef association_matrix< BlockType, row_grouped< 1 > > space_type;
    typedef typename space_type::block_type                 block_type;
    typedef size_t *                                        result_type;

    typedef typename space_type::raw_block_pointer    block_pointer;

    typedef clotho::utility::debruijn_bit_walker< block_type >  bit_walker_type;
    typedef clotho::utility::BitHelper< block_type >            bit_helper_type;

    void operator()( space_type & ss, result_type res, size_t & fixed_offset, size_t & lost_offset, size_t & free_count, size_t M ) {

        size_t W = bit_helper_type::padded_block_count( M );
#ifdef  DEBUGGING
        std::cerr << "Free Space Padded size: " << W << std::endl;
#endif  // DEBUGGING

        resize( W );

        memset( tmp, 255, sizeof(block_type) * W );
        memset( tmp + W, 0, sizeof(block_type) * W );

        const size_t row_count = ss.block_row_count();

        size_t fo = fixed_offset;
        size_t lo = lost_offset;
        size_t fr = free_count;

        size_t k = 0;
        while( k < row_count ) {

            block_pointer start = ss.begin_block_row( k );
            unsigned int N = ss.getSoftSize(k);
            if( N ) {
//                updateFixed( start, tmp, N );
//                updateVariable( start, tmp + W , N );
                updateFixedVariable( start, tmp, tmp + W, N );
            }

            for( ; N < W; ++N )
                tmp[ N ] = bit_helper_type::ALL_UNSET;
            ++k;
        }

        size_t j = 0;
        
        for( unsigned int i = 0; i < W; ++i ) {
            block_type fx = tmp[ i ];
            block_type var = tmp[ W + i ];

            block_type ls = ~(fx | var);

            while( fx ) {
                size_t b_idx = bit_walker_type::unset_next_index( fx ) + j;
                if( b_idx < M ) {
                    res[ fr++ ] = b_idx;
                    res[ fo++ ] = b_idx;
                }
            }

            while( ls ) {
                size_t idx = bit_walker_type::unset_next_index( ls ) + j;
                if( idx < M ) {
                    res[ fr++ ] = idx;
                    res[ lo++ ] = idx;
                }
            }

            j += bit_helper_type::BITS_PER_BLOCK;
        }

        fixed_offset = fo;
        lost_offset = lo;
        free_count = fr;
    }

    virtual ~free_space_evaluator() {
        if( tmp != NULL ) {
            delete [] tmp;
        }
    }

protected:

    void updateFixedVariable( block_type * start, block_type * fx_out, block_type * var_out, unsigned int N ) {
        assert( (N % bit_helper_type::BLOCKS_PER_SSE) == 0 );

        __m128i * s = reinterpret_cast< __m128i * >( start );
        __m128i * o = reinterpret_cast< __m128i * >( fx_out );
        __m128i * vo = reinterpret_cast< __m128i * >( var_out );

        unsigned int M = N / bit_helper_type::BLOCKS_PER_SSE;
        while( M-- ) {
            __m128i b = *s++;
            *o = (*o & b);
            ++o;
            *vo = (*vo | b);
            ++vo;
        }
    }

    inline void updateFixed( block_type * start, block_type * out, unsigned int N ) {
        assert( (N % bit_helper_type::BLOCKS_PER_SSE) == 0 );
        unsigned int i = 0;
        __m128i * s = reinterpret_cast< __m128i * >( start );
        __m128i * o = reinterpret_cast< __m128i * >( out );
        unsigned int M = N / bit_helper_type::BLOCKS_PER_SSE;
        while( i < M ) {
            o[ i ] = (o[ i ] & s[ i ]);
            ++i;
        }
    }

    inline void updateVariable( block_type * start, block_type * out, unsigned int N ) {
        assert( (N % bit_helper_type::BLOCKS_PER_SSE) == 0 );
        unsigned int i = 0;
        __m128i * s = reinterpret_cast< __m128i * >( start );
        __m128i * o = reinterpret_cast< __m128i * >( out );
        unsigned int M = N / bit_helper_type::BLOCKS_PER_SSE;
        while( i < M ) {
            o[ i ] = (o[ i ] | s[i]);
            ++i;
        }
    }

    void resize( size_t W ) {
        if( 2 * W > m_alloc_size ) {
            if( tmp != NULL ) {
                delete [] tmp;
            }

            tmp = new block_type[ 2 * W ];
            m_alloc_size = 2 * W;
        }
    }
    block_type * tmp;
    size_t m_alloc_size;
};
}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_ROW_GROUPED_FREE_SPACE_EVALUATOR_HPP_
