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

        size_t k = 0;
        while( k < row_count ) {

            block_pointer start = ss.begin_block_row( k );
            block_pointer end = ss.end_block_row( k );

            size_t fix_idx = 0, var_idx = W;
            while( start != end ) {
                block_type fx = bit_helper_type::ALL_SET;
                block_type var = bit_helper_type::ALL_UNSET;

                size_t i = row_grouped< P >::GROUP_SIZE;
                while( i-- ) {
                    block_type b = *start++;
                    fx &= b;
                    var |= b;
                }
                tmp[ fix_idx++ ] &= fx;
                tmp[ var_idx++ ] |= var;
            }

            k += row_grouped< P >::GROUP_SIZE;
        }

        size_t fix_idx = 0, var_idx = W, j = 0;
        while( fix_idx < W ) {
            block_type fx = tmp[ fix_idx++ ];
            block_type var = tmp[ var_idx++ ];

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

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_ROW_GROUPED_FREE_SPACE_EVALUATOR_HPP_
