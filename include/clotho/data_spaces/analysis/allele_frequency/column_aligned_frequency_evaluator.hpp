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
#ifndef CLOTHO_COLUMN_ALIGNED_FREQUENCY_EVALUATOR_HPP_
#define CLOTHO_COLUMN_ALIGNED_FREQUENCY_EVALUATOR_HPP_

#include "clotho/data_spaces/analysis/allele_frequency/frequency_evaluator_def.hpp"
#include "clotho/data_spaces/association_matrix/column_aligned_association_matrix.hpp"
#include <boost/dynamic_bitset.hpp>
#include "clotho/utility/debruijn_bit_walker.hpp"

namespace clotho {
namespace genetics {

template < class BlockType >
struct frequency_evaluator< association_matrix< BlockType, column_aligned > > {
    typedef association_matrix< BlockType, column_aligned > space_type;
    typedef typename space_type::block_type                 block_type;

    typedef typename space_type::bit_helper_type                        bit_helper_type;
    typedef typename clotho::utility::debruijn_bit_walker< block_type > bit_walker_type;

    typedef size_t *                                        result_type;

    void operator()( space_type & ss, boost::dynamic_bitset<> & indices, result_type * res ) {
        size_t M = ss.column_count();
        size_t N = ss.block_column_count();

        size_t buffer[ bit_helper_type::BITS_PER_BLOCK ];
        memset( buffer, 0, bit_helper_type::BITS_PER_BLOCK * sizeof( size_t ) );

        typedef typename space_type::raw_block_pointer iterator;

        size_t i = 0, j = M;
        while( j ) {
            iterator first = ss.begin_block_row( i );
            iterator last = ss.end_block_row( i );

            size_t k = 0;
            while( first != last ) {
                block_type b = *first++;
                if( indices.test(k++) ) {
                    while( b ) {
                        unsigned int b_idx = bit_walker_type::unset_next_index( b );
                        buffer[ b_idx ] += 1;
                    }
                }
            }

            size_t l = (( j < bit_helper_type::BITS_PER_BLOCK) ? j : bit_helper_type::BITS_PER_BLOCK);

            memcpy( res + i * bit_helper_type::BITS_PER_BLOCK, buffer, sizeof( size_t ) * l );
            memset( buffer, 0, bit_helper_type::BITS_PER_BLOCK * sizeof( size_t ) );

            j -= l;
            ++i;
        }
    } 
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_COLUMN_ALIGNED_FREQUENCY_EVALUATOR_HPP_
