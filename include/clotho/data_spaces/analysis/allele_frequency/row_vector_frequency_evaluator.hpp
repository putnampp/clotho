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
#ifndef CLOTHO_ROW_VECTOR_FREQUENCY_EVALUATOR_HPP_
#define CLOTHO_ROW_VECTOR_FREQUENCY_EVALUATOR_HPP_

#include "clotho/data_spaces/analysis/allele_frequency/frequency_evaluator_def.hpp"
#include "clotho/data_spaces/association_matrix/row_vector_association_matrix.hpp"
#include <boost/dynamic_bitset.hpp>
#include "clotho/utility/debruijn_bit_walker.hpp"

namespace clotho {
namespace genetics {

template < class BlockType >
struct frequency_evaluator< association_matrix< BlockType, row_vector > > {
    typedef association_matrix< BlockType, row_vector > space_type;
    typedef typename space_type::block_type                 block_type;
    typedef size_t *                                        result_type;

    typedef typename space_type::bit_helper_type                        bit_helper_type;
    typedef typename clotho::utility::debruijn_bit_walker< block_type > bit_walker_type;

    void operator()( space_type & ss, boost::dynamic_bitset<> & indices, result_type res ) {
        typedef typename space_type::raw_block_pointer iterator;

        size_t N = ss.row_count();
        size_t i = 0;
        while ( i < N ) {
            if( !indices.test(i) ) {
                ++i;
                continue;
            }

            iterator start = ss.begin_row(i);
            iterator end = ss.end_row(i);

            size_t j = 0;
            while( start != end ) {
                block_type b = *start++;

                while( b ) {
                    unsigned int b_idx = bit_walker_type::unset_next_index( b ) + j;
                    res[ b_idx ] += 1;
                }
                
                j += bit_helper_type::BITS_PER_BLOCK;
            }
            ++i;
        }
    }
};

}   // namespace genetics
}   // namespace clotho
#endif  // CLOTHO_ROW_VECTOR_FREQUENCY_EVALUATOR_HPP_

