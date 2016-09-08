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
#ifndef CLOTHO_ROW_VECTOR_FREE_SPACE_EVALUATOR_HPP_
#define CLOTHO_ROW_VECTOR_FREE_SPACE_EVALUATOR_HPP_

#include "clotho/data_spaces/association_matrix/row_vector_association_matrix.hpp"
#include "clotho/utility/bit_helper.hpp"
#include "clotho/utility/debruijn_bit_walker.hpp"
#include <map>
#include <set>

namespace clotho {
namespace genetics {


template < class BlockType, class SizeType >
class free_space_evaluator< association_matrix< BlockType, row_vector >, SizeType > {
public:
    typedef association_matrix< BlockType, row_vector >     space_type;
    typedef typename space_type::raw_vector                 raw_vector;
    typedef typename raw_vector::block_type                 block_type;

    typedef SizeType                                        size_type;
    typedef size_type *                                        result_type;

    typedef typename raw_vector::raw_pointer                raw_pointer;

    typedef clotho::utility::debruijn_bit_walker< block_type >  bit_walker_type;
    typedef clotho::utility::BitHelper< block_type >            bit_helper_type;

    free_space_evaluator() :
        tmp(NULL)
        , m_alloc_size(0)
    {}

    void operator()( space_type & ss, result_type res, size_type & fixed_offset, size_type & lost_offset, size_type & free_count, size_type M ) {

        size_type W = ss.hard_block_count();

        resize( W );

        for( size_type i = 0; i < W; ++i ) {
            tmp[ 2 * i ] = bit_helper_type::ALL_SET;
            tmp[ 2 * i + 1 ] = bit_helper_type::ALL_UNSET;
        }

        std::set< raw_pointer > analyzed;

        const size_type row_count = ss.row_count();

        size_type fo = fixed_offset;
        size_type lo = lost_offset;
        size_type fr = free_count;
        
        size_type k = 0;
        while( k < row_count ) {

            raw_pointer start = ss.getRow( k ).get();

            if( analyzed.find( start ) == analyzed.end() ) {
                analyzed.insert( start );
                size_type  N = ss.getRow( k ).m_size;

#ifdef DEBUGGING
                assert( N <= W );
#endif // DEBUGGING

                raw_pointer end = start + N;
                block_type * t = tmp;
                while( start != end ) {
                    block_type b = *start++;
                    *t++ &= b;
                    *t++ |= b;
                }

                size_type i = N;
                while( i < W ) {
                    tmp[ 2 * i ] = bit_helper_type::ALL_UNSET;
                    ++i;
                }
            }

            ++k;
        }

#ifdef DEBUGGING
        std::cerr << "Free Space analyzed: " << analyzed.size() << std::endl;
#endif  // DEBUGGING

        size_type j = 0;
        
        for( unsigned int i = 0; i < W; ++i ) {
            block_type fx = tmp[ 2 *i ];
            block_type var = tmp[2 * i + 1];

            block_type ls = ~(fx | var);

            while( fx ) {
                size_type b_idx = bit_walker_type::unset_next_index( fx ) + j;
                if( b_idx < M ) {
                    res[ fr++ ] = b_idx;
                    res[ fo++ ] = b_idx;
                }
            }

            while( ls ) {
                size_type idx = bit_walker_type::unset_next_index( ls ) + j;
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
    void resize( size_type W ) {
        if( 2 * W > m_alloc_size ) {
            if( tmp != NULL ) {
                delete [] tmp;
            }

            tmp = new block_type[ 2 * W ];
            m_alloc_size = 2 * W;
        }
    }
    block_type * tmp;
    size_type m_alloc_size;
};

/*
/// this approach eliminates writing to the heap
/// utilizes stack buffers to write temporary data
template < class BlockType >
class free_space_evaluator< association_matrix< BlockType, row_vector > > {
public:
    typedef association_matrix< BlockType, row_vector > space_type;
    typedef typename space_type::raw_vector                 raw_vector;
    typedef typename raw_vector::block_type                 block_type;
    typedef size_type *                                        result_type;

    typedef typename raw_vector::raw_pointer                raw_pointer;

    typedef clotho::utility::debruijn_bit_walker< block_type >  bit_walker_type;
    typedef clotho::utility::BitHelper< block_type >            bit_helper_type;

    free_space_evaluator() {}

    void operator()( space_type & ss, result_type res, size_type & fixed_offset, size_type & lost_offset, size_type & free_count, size_type M ) {

        size_type W = ss.hard_block_count();

        std::map< raw_pointer, size_type > analyzed;
        typedef typename std::map< raw_pointer, size_type >::iterator iterator;

        const size_type row_count = ss.row_count();

        size_type fo = fixed_offset;
        size_type lo = lost_offset;
        size_type fr = free_count;
        
        // pre-scan population for all unique sequences
        size_type k = 0;
        while( k < row_count ) {
            raw_pointer start = ss.getRow( k ).get();
            size_type s = ss.getRow( k ).m_size;

            if( analyzed.find( start ) == analyzed.end() ) {
                analyzed[ start ] = s;
            }
            ++k;
        }

#ifdef DEBUGGING
        std::cerr << "Free Space analyzing: " << analyzed.size() << std::endl;
#endif  // DEBUGGING

        const size_type BUFFER_SIZE = 128;
        k = 0;
        while( k < W ) {
//            // 16 * sizeof(block_type) * 2 == 16 * 8 * 2 == 256 byte buffer
            block_type fx_buffer[ BUFFER_SIZE ];
            block_type var_buffer[ BUFFER_SIZE ];

            // initialize buffers
            for( int x = 0; x < BUFFER_SIZE; ++x ) {
                fx_buffer[ x ] = bit_helper_type::ALL_SET;
                var_buffer[ x ] = bit_helper_type::ALL_UNSET;
            }

            // analyze block columns
            for( iterator it = analyzed.begin(); it != analyzed.end(); it++ ) {
                size_type S = it->second;
                size_type N = ((k >= S)? 0 : (( k + BUFFER_SIZE <= S ) ? BUFFER_SIZE : (S - k)));
                size_type x = 0;
                raw_pointer r = it->first;
                while( x < N ) {
                    block_type b = r[ k + x ];
                    fx_buffer[ x ] &= b;
                    var_buffer[ x ] |= b;
                    ++x;
                }

                while( x < BUFFER_SIZE ) {
                    fx_buffer[ x ] = bit_helper_type::ALL_UNSET;
                    ++x;
                }
            }


            // write results
            size_type j = k * bit_helper_type::BITS_PER_BLOCK;

            for( unsigned int i = 0; i < BUFFER_SIZE; ++i ) {
                block_type fx = fx_buffer[ i ];
                block_type var = var_buffer[ i ];

                block_type ls = ~(fx | var);

                while( fx ) {
                    size_type b_idx = bit_walker_type::unset_next_index( fx ) + j;
                    if( b_idx < M ) {
                        res[ fr++ ] = b_idx;
                        res[ fo++ ] = b_idx;
                    }
                }

                while( ls ) {
                    size_type idx = bit_walker_type::unset_next_index( ls ) + j;
                    if( idx < M ) {
                        res[ fr++ ] = idx;
                        res[ lo++ ] = idx;
                    }
                }

                j += bit_helper_type::BITS_PER_BLOCK;
            }
            
            k += BUFFER_SIZE;
        }

        fixed_offset = fo;
        lost_offset = lo;
        free_count = fr;
    }

    virtual ~free_space_evaluator() { }
};
*/
}   // namespace genetics
}   // namespace clotho
#endif  // CLOTHO_ROW_VECTOR_FREE_SPACE_EVALUATOR_HPP_
