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
#ifndef CLOTHO_POPULATION_SPACE_FREQUENCY_EVALUATOR_HPP_
#define CLOTHO_POPULATION_SPACE_FREQUENCY_EVALUATOR_HPP_

#include "clotho/data_spaces/analysis/allele_frequency/frequency_evaluator_def.hpp"
#include "clotho/data_spaces/population_space/population_space.hpp"
#include <boost/dynamic_bitset.hpp>
#include "clotho/utility/debruijn_bit_walker.hpp"

namespace clotho {
namespace genetics {

template < class BlockType, class WeightType, class SizeType >
struct frequency_evaluator< population_space< BlockType, WeightType >, SizeType > {
    typedef population_space< BlockType, WeightType >   space_type;
    typedef typename space_type::base_genome_type::sequence_type::block_type block_type;

    typedef SizeType                                        size_type;
    typedef size_type *                                     result_type;

    typedef typename space_type::haploid_genomes            haploid_genomes;

    typedef typename space_type::base_genome_type::sequence_type::bit_helper_type   bit_helper_type;
    typedef typename clotho::utility::debruijn_bit_walker< block_type >             bit_walker_type;

    void operator()( const space_type & ss, const boost::dynamic_bitset<> & indices, result_type column_margin, result_type row_margin ) {
        haploid_genomes local_genomes;
        buildGenomes( ss, indices, local_genomes );

        unsigned int i = 0;
        for( typename haploid_genomes::const_iterator it = local_genomes.begin(); it != local_genomes.end(); ++it ) {
            if( it->first ) {
                const size_type N = it->second;

                unsigned int M = 0, j = 0;
                typename space_type::const_sequence_iterator first = it->first->begin_sequence(), last = it->first->end_sequence();
                while( first != last ) {
                    block_type b = *first++;
                    while( b ) {
                        unsigned int b_idx = bit_walker_type::unset_next_index( b ) + j;
                        column_margin[ b_idx ] += N;
                        ++M;
                    }
                    j += bit_helper_type::BITS_PER_BLOCK;
                }
                row_margin[ i++ ] = M;
            }
        }
    }

    void buildGenomes( const space_type & ss, const boost::dynamic_bitset<> & indices, haploid_genomes & local_genomes ) {
        unsigned int idx = 0;
        for( typename space_type::const_individual_iterator it = ss.begin_individual(); it != ss.end_individual(); ++it ) {
            if( indices[ idx++ ] ) {
                if( local_genomes.find( it->first) == local_genomes.end() ) {
                    local_genomes.insert( std::make_pair( it->first, 1 ) );
                } else {
                    local_genomes[ it->first ] += 1;
                }
            }

            if( indices[ idx++ ] ) {
                if( local_genomes.find( it->second) == local_genomes.end() ) {
                    local_genomes.insert( std::make_pair( it->second, 1 ) );
                } else {
                    local_genomes[ it->second ] += 1;
                }

            }
        }

    }
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_POPULATION_SPACE_FREQUENCY_EVALUATOR_HPP_
