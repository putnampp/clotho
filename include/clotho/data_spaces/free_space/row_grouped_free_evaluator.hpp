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
#include <map>
#include <set>

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
/*
template < class BlockType >
class free_space_evaluator< association_matrix< BlockType, row_grouped< 1 > > > {
public:
    typedef association_matrix< BlockType, row_grouped< 1 > > space_type;
    typedef typename space_type::row_vector                 row_vector;
    typedef typename row_vector::block_type                 block_type;
    typedef size_t *                                        result_type;

    typedef typename row_vector::raw_pointer                raw_pointer;

    typedef clotho::utility::debruijn_bit_walker< block_type >  bit_walker_type;
    typedef clotho::utility::BitHelper< block_type >            bit_helper_type;

    free_space_evaluator() :
        tmp(NULL)
        , m_alloc_size(0)
    {}

    void operator()( space_type & ss, result_type res, size_t & fixed_offset, size_t & lost_offset, size_t & free_count, size_t M ) {

//        size_t W = bit_helper_type::padded_block_count( M );
        size_t W = ss.hard_block_count();

        resize( W );

        for( size_t i = 0; i < W; ++i ) {
            tmp[ 2 * i ] = bit_helper_type::ALL_SET;
            tmp[ 2 * i + 1 ] = bit_helper_type::ALL_UNSET;
        }

        std::set< raw_pointer > analyzed;

        const size_t row_count = ss.row_count();

        size_t fo = fixed_offset;
        size_t lo = lost_offset;
        size_t fr = free_count;
        
        size_t k = 0;
        while( k < row_count ) {

            raw_pointer start = ss.getRow( k ).get();

            if( analyzed.find( start ) == analyzed.end() ) {
                analyzed.insert( start );
                size_t  N = ss.getRow( k ).m_size;

#ifdef DEBUGGING
                assert( N <= W );
#endif // DEBUGGING

//                size_t i = 0;
//                while( i < N ) {
//                    block_type b = start[i];
//                    tmp[ 2 * i ] &= b;
//                    tmp[ 2 * i + 1 ] |= b;
//                    ++i;
//                }

                raw_pointer end = start + N;
                block_type * t = tmp;
                while( start != end ) {
                    *t++ &= *start;
                    *t++ |= *start++;
                }

                size_t i = N;
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

        size_t j = 0;
        
        for( unsigned int i = 0; i < W; ++i ) {
            block_type fx = tmp[ 2 *i ];
            block_type var = tmp[2 * i + 1];

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
};*/

/// this approach eliminates writing to the heap
/// utilizes stack buffers to write temporary data
template < class BlockType >
class free_space_evaluator< association_matrix< BlockType, row_grouped< 1 > > > {
public:
    typedef association_matrix< BlockType, row_grouped< 1 > > space_type;
    typedef typename space_type::row_vector                 row_vector;
    typedef typename row_vector::block_type                 block_type;
    typedef size_t *                                        result_type;

    typedef typename row_vector::raw_pointer                raw_pointer;

    typedef clotho::utility::debruijn_bit_walker< block_type >  bit_walker_type;
    typedef clotho::utility::BitHelper< block_type >            bit_helper_type;

    free_space_evaluator() {}

    void operator()( space_type & ss, result_type res, size_t & fixed_offset, size_t & lost_offset, size_t & free_count, size_t M ) {

        size_t W = ss.hard_block_count();

        std::map< raw_pointer, size_t > analyzed;
        typedef typename std::map< raw_pointer, size_t >::iterator iterator;

        const size_t row_count = ss.row_count();

        size_t fo = fixed_offset;
        size_t lo = lost_offset;
        size_t fr = free_count;
        
        // pre-scan population for all unique sequences
        size_t k = 0;
        while( k < row_count ) {
            raw_pointer start = ss.getRow( k ).get();
            size_t s = ss.getRow( k ).m_size;

            if( analyzed.find( start ) == analyzed.end() ) {
                analyzed[ start ] = s;
            }
            ++k;
        }

#ifdef DEBUGGING
        std::cerr << "Free Space analyzing: " << analyzed.size() << std::endl;
#endif  // DEBUGGING

        const size_t BUFFER_SIZE = 8;
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
                size_t S = it->second;
                size_t N = ((k >= S)? 0 : (( k + BUFFER_SIZE <= S ) ? BUFFER_SIZE : (S - k)));
                size_t x = 0;
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
            size_t j = k * bit_helper_type::BITS_PER_BLOCK;

            for( unsigned int i = 0; i < BUFFER_SIZE; ++i ) {
                block_type fx = fx_buffer[ i ];
                block_type var = var_buffer[ i ];

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
            
            k += BUFFER_SIZE;
        }

        fixed_offset = fo;
        lost_offset = lo;
        free_count = fr;
    }

    virtual ~free_space_evaluator() { }
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_ROW_GROUPED_FREE_SPACE_EVALUATOR_HPP_
