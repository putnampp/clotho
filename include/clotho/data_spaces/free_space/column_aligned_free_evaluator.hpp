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
#ifndef CLOTHO_COLUMN_ALIGNED_FREE_SPACE_EVALUATOR_HPP_
#define CLOTHO_COLUMN_ALIGNED_FREE_SPACE_EVALUATOR_HPP_

#include "clotho/data_spaces/association_matrix/column_aligned_association_matrix.hpp"
#include "clotho/utility/bit_helper.hpp"
#include "clotho/utility/debruijn_bit_walker.hpp"

namespace clotho {
namespace genetics {

template < class BlockType >
class free_space_evaluator< association_matrix< BlockType, column_aligned > > {
public:
    typedef association_matrix< BlockType, column_aligned > space_type;
    typedef typename space_type::block_type                 block_type;
    typedef size_t *                                        result_type;

    typedef typename space_type::raw_block_pointer    block_pointer;

    typedef clotho::utility::debruijn_bit_walker< block_type >  bit_walker_type;

    typedef clotho::utility::BitHelper< block_type >    bit_helper_type;

    void operator()( space_type & ss, result_type res, size_t & fixed_offset, size_t & lost_offset, size_t & free_count, size_t M ) {
        const size_t col_count = ss.block_column_count();
        const size_t row_count = ss.block_row_count();

        const size_t N = col_count % 4;

        size_t j = 0, k = 0;

        size_t fo = fixed_offset;
        size_t lo = lost_offset;
        size_t fr = free_count;

        while( k < row_count ) {
            block_type fx = bit_helper_type::ALL_SET, var = bit_helper_type::ALL_UNSET;

            block_pointer start = ss.begin_block_row( k );
            block_pointer end = ss.end_block_row( k );

            // unwinding loop
            size_t i = N;
            while( i-- ) {
                block_type v = *start++;
                fx &= v;
                var |= v;
            }

            while( start != end ) {
                block_type v = *start++;
                fx &= v;
                var |= v;

                v = *start++;
                fx &= v;
                var |= v;

                v = *start++;
                fx &= v;
                var |= v;

                v = *start++;
                fx &= v;
                var |= v;
            }

            block_type ls = ~(fx | var);

            while( fx ) {
                size_t idx = bit_walker_type::unset_next_index( fx ) + j;
                if( idx < M ) {
                    res[ fr++ ] = idx;
                    res[ fo++ ] =  idx;
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
            ++k;
        }

        fixed_offset = fo;
        lost_offset = lo;
        free_count = fr;
    }
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_COLUMN_ALIGNED_FREE_SPACE_EVALUATOR_HPP_
