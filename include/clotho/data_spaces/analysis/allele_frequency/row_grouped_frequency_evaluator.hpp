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
#ifndef CLOTHO_ROW_GROUPED_FREQUENCY_EVALUATOR_HPP_
#define CLOTHO_ROW_GROUPED_FREQUENCY_EVALUATOR_HPP_

#include "clotho/data_spaces/analysis/allele_frequency/frequency_evaluator_def.hpp"
#include "clotho/data_spaces/association_matrix/row_grouped_association_matrix.hpp"
#include <boost/dynamic_bitset.hpp>
#include "clotho/utility/debruijn_bit_walker.hpp"

namespace clotho {
namespace genetics {

template < class BlockType, class SizeType, unsigned char P >
struct frequency_evaluator< association_matrix< BlockType, row_grouped< P > >, SizeType > {
    typedef association_matrix< BlockType, row_grouped< P > > space_type;
    typedef typename space_type::block_type                 block_type;
    typedef SizeType                                        size_type;
    typedef size_type *                                        result_type;

    typedef typename space_type::bit_helper_type                        bit_helper_type;
    typedef typename clotho::utility::debruijn_bit_walker< block_type > bit_walker_type;

    void operator()( space_type & ss, boost::dynamic_bitset<> & indices, result_type column_margin, result_type row_margin ) {
        typedef typename space_type::raw_block_pointer iterator;

        size_type N = ss.row_count();
        size_type i = 0;
        while ( i < N ) {
            if( !indices.test(i) ) {
                ++i;
                continue;
            }

            iterator start = ss.begin_row(i);
            iterator end = ss.end_row(i);

            size_type j = 0, M = 0;
            while( start != end ) {
                block_type b = *start;

                while( b ) {
                    unsigned int b_idx = bit_walker_type::unset_next_index( b ) + j;
                    column_margin[ b_idx ] += 1;
                    ++M;
                }
                
                start += row_grouped< P >::GROUP_SIZE;
                j += bit_helper_type::BITS_PER_BLOCK;
            }
            row_margin[ i ] = M;
            ++i;
        }
    }
};

template < class BlockType, class SizeType >
struct frequency_evaluator< association_matrix< BlockType, row_grouped< 1 > >, SizeType > {
    typedef association_matrix< BlockType, row_grouped< 1 > > space_type;
    typedef typename space_type::block_type                 block_type;
    typedef SizeType                                        size_type;
    typedef size_type *                                     result_type;

    typedef typename space_type::bit_helper_type                        bit_helper_type;
    typedef typename clotho::utility::debruijn_bit_walker< block_type > bit_walker_type;

    void operator()( space_type & ss, boost::dynamic_bitset<> & indices, result_type column_margin, result_type row_margin ) {
        typedef typename space_type::raw_block_pointer iterator;

        size_type N = ss.row_count();
        size_type i = 0;
        while ( i < N ) {
            if( !indices.test(i) ) {
                ++i;
                continue;
            }

            iterator start = ss.begin_row(i);
            iterator end = ss.end_row(i);

            size_type j = 0, M = 0;
            while( start != end ) {
                block_type b = *start++;

                while( b ) {
                    unsigned int b_idx = bit_walker_type::unset_next_index( b ) + j;
                    column_margin[ b_idx ] += 1;
                    ++M;
                }
                
                j += bit_helper_type::BITS_PER_BLOCK;
            }
            row_margin[ i ] = M;
            ++i;
        }
    }
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_ROW_GROUPED_FREQUENCY_EVALUATOR_HPP_
