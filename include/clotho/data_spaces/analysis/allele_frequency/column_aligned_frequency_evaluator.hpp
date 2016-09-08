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
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either expcolumn_margins or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License.
#ifndef CLOTHO_COLUMN_ALIGNED_FREQUENCY_EVALUATOR_HPP_
#define CLOTHO_COLUMN_ALIGNED_FREQUENCY_EVALUATOR_HPP_

#include "clotho/data_spaces/analysis/allele_frequency/frequency_evaluator_def.hpp"
#include "clotho/data_spaces/association_matrix/column_aligned_association_matrix.hpp"
#include <boost/dynamic_bitset.hpp>
#include "clotho/utility/debruijn_bit_walker.hpp"

namespace clotho {
namespace genetics {

template < class BlockType, class SizeType >
struct frequency_evaluator< association_matrix< BlockType, column_aligned >, SizeType > {
    typedef association_matrix< BlockType, column_aligned > space_type;
    typedef typename space_type::block_type                 block_type;

    typedef typename space_type::bit_helper_type                        bit_helper_type;
    typedef typename clotho::utility::debruijn_bit_walker< block_type > bit_walker_type;

    typedef SizeType                                        size_type;
    typedef size_type *                                        result_type;

    void operator()( space_type & ss, boost::dynamic_bitset<> & indices, result_type column_margin, result_type row_margin ) {
        size_type M = ss.column_count();
//        size_type N = ss.block_column_count();

        size_type buffer[ bit_helper_type::BITS_PER_BLOCK ];
        memset( buffer, 0, bit_helper_type::BITS_PER_BLOCK * sizeof( size_type ) );

        typedef typename space_type::raw_block_pointer iterator;

        size_type i = 0, j = M;
        while( j ) {
            iterator first = ss.begin_block_row( i );
            iterator last = ss.end_block_row( i );

            size_type k = 0;
            while( first != last ) {
                block_type b = *first++;
                if( indices.test(k) ) {
                    size_type m = 0;
                    while( b ) {
                        unsigned int b_idx = bit_walker_type::unset_next_index( b );
                        buffer[ b_idx ] += 1;
                        ++m;
                    }
                    row_margin[ k ] += m;
                }
                ++k;
            }

            size_type l = (( j < bit_helper_type::BITS_PER_BLOCK) ? j : bit_helper_type::BITS_PER_BLOCK);

            memcpy( column_margin + i * bit_helper_type::BITS_PER_BLOCK, buffer, sizeof( size_type ) * l );
            memset( buffer, 0, bit_helper_type::BITS_PER_BLOCK * sizeof( size_type ) );

            j -= l;
            ++i;
        }
    } 
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_COLUMN_ALIGNED_FREQUENCY_EVALUATOR_HPP_
