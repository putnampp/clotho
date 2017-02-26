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
#ifndef CLOTHO_FREQUENCY_EVALUATOR_POPULATION_SPACE_ROW_ALIGNMENT_HPP_
#define CLOTHO_FREQUENCY_EVALUATOR_POPULATION_SPACE_ROW_ALIGNMENT_HPP_

#include "clotho/data_spaces/analysis/allele_frequency/frequency_evaluator_def.hpp"
#include "clotho/data_spaces/population_space/population_spaces.hpp"
#include <boost/dynamic_bitset.hpp>
#include "clotho/utility/debruijn_bit_walker.hpp"

namespace clotho {
namespace genetics {

template < class BlockType, class WeightType, class SizeType >
struct frequency_evaluator< population_space< row_block_alignment< BlockType >, trait_space_vector< WeightType > >, SizeType > {
    typedef population_space< row_block_alignment< BlockType >, trait_space_vector< WeightType > >   space_type;
    typedef typename space_type::block_type                 block_type;
    typedef typename space_type::genome_pointer             genome_pointer;

    typedef SizeType                                        size_type;
    typedef size_type *                                     result_type;

    typedef typename space_type::bit_helper_type            bit_helper_type;
    typedef typename clotho::utility::debruijn_bit_walker< block_type >             bit_walker_type;

    void operator()( space_type & ss, boost::dynamic_bitset<> & indices, result_type column_margin, result_type row_margin ) {

        const unsigned int N = ss.haploid_genome_count();

        for( unsigned int i = 0; i < N; ++i ) {
            if( !indices.test(i) ) continue;

            genome_pointer first = ss.begin_genome(i), last = ss.end_genome(i);

            unsigned int j = 0;
            unsigned int M = 0;
            while( first != last ) {
                block_type b = *first++;

                unsigned int idx = j;
                while( b ) {
                    idx += bit_walker_type::next_and_shift( b );
                    column_margin[ idx ]++;
                    ++M;
                }
                j += bit_helper_type::BITS_PER_BLOCK;
            }

            row_margin[ i ] = M;
        }
    }

    virtual ~frequency_evaluator() {}

};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_FREQUENCY_EVALUATOR_POPULATION_SPACE_ROW_ALIGNMENT_HPP_
