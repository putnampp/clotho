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
#ifndef CLOTHO_POPULATION_SPACE_COLUMNAR_FREQUENCY_EVALUATOR_HPP_
#define CLOTHO_POPULATION_SPACE_COLUMNAR_FREQUENCY_EVALUATOR_HPP_

#include "clotho/data_spaces/analysis/allele_frequency/frequency_evaluator_def.hpp"
#include "clotho/data_spaces/population_space/population_space_columnar.hpp"
#include <boost/dynamic_bitset.hpp>
#include "clotho/utility/debruijn_bit_walker.hpp"

namespace clotho {
namespace genetics {

template < class BlockType, class WeightType, class SizeType >
struct frequency_evaluator< population_space_columnar< BlockType, WeightType >, SizeType > {
    typedef population_space_columnar< BlockType, WeightType >   space_type;
    typedef typename space_type::block_type block_type;

    typedef SizeType                                        size_type;
    typedef size_type *                                     result_type;

    typedef typename space_type::bit_helper_type            bit_helper_type;
    typedef typename clotho::utility::debruijn_bit_walker< block_type >             bit_walker_type;

    void operator()( space_type & ss, boost::dynamic_bitset<> & indices, result_type column_margin, result_type row_margin ) {

        for( unsigned int i = 0; i < ss.getMaxBlocks(); ++i ) {
            typename space_type::row_pointer first = ss.begin_block_row( i ), last = ss.end_block_row( i );

            const unsigned int j = i * bit_helper_type::BITS_PER_BLOCK;
            unsigned int k = 0;
            while( first != last ) {
                if( indices.test( k ) ) {
                    block_type b = *first;

                    unsigned int idx = j;
                    unsigned int M = 0;
                    while( b ) {
                        idx += bit_walker_type::next_and_shift( b );
                        column_margin[ idx ]++;
                        ++M;
                    }
                    row_margin[ k ] += M;
                }

                ++k;
                ++first;
            }
        }
    }

    virtual ~frequency_evaluator() {}

};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_POPULATION_SPACE_COLUMNAR_FREQUENCY_EVALUATOR_HPP_
