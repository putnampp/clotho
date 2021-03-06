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
#ifndef CLOTHO_ROW_VECTOR_WEIGHT_ACCUMULATOR_HPP_
#define CLOTHO_ROW_VECTOR_WEIGHT_ACCUMULATOR_HPP_

#include "clotho/data_spaces/association_matrix/row_vector_association_matrix.hpp"
#include <iostream>

#include "clotho/utility/debruijn_bit_walker.hpp"

#include <unordered_map>

namespace clotho {
namespace genetics {

template < class AlleleSpace, class BlockType >
class weight_accumulator< AlleleSpace, association_matrix< BlockType, row_vector > > {
public:
    typedef AlleleSpace                                         allele_type;
    typedef association_matrix< BlockType, row_vector >     association_type;

    typedef typename allele_type::weight_type                   weight_type;
    typedef typename allele_type::weight_pointer                weight_pointer;

    typedef typename association_type::raw_vector               raw_vector;
    typedef typename raw_vector::block_type                     block_type;
    typedef typename association_type::raw_block_pointer        block_pointer;

    typedef clotho::utility::debruijn_bit_walker< block_type >  bit_walker_type;
    

    void operator()( allele_type & alls, association_type & seqs, weight_pointer res, size_t trait_count ) {
        size_t  M = seqs.row_count();

        std::unordered_map< block_pointer, size_t > analyzed;
        
        size_t i = 0;
        while( i < M ) {

            block_pointer start = seqs.begin_row( i );

            if( analyzed.find( start ) == analyzed.end() ) {
                analyzed[ start ] = i;

                block_pointer end = seqs.end_row(i);

                size_t j = 0;
                while( start != end ) {
                    block_type b = *start++;

                    while( b ) {
                        size_t b_idx = j + bit_walker_type::unset_next_index( b );

                        weight_pointer t_start = alls.begin_trait_weight( b_idx );
                        weight_pointer t_end = alls.end_trait_weight( b_idx );

                        size_t a = i * trait_count;
                        while( t_start != t_end ) {
                            res[ a++ ] += *t_start++;
                        }
                    }
                    
                    j += association_type::bit_helper_type::BITS_PER_BLOCK;
                }

            } else {
                size_t src = analyzed[ start ] * trait_count;
                size_t a = i * trait_count;

                size_t x = 0;
                while( x < trait_count ) {
                    res[ a + x ] = res[ src + x ];
                    ++x;
                }
            }
            ++i;
        }

#ifdef DEBUGGING
        std::cerr << "Weights accumulated for: " << analyzed.size() << std::endl;
#endif  // DEBUGGING
    }
};
}   // namespace genetics
}   // namespace clotho
#endif  // CLOTHO_ROW_VECTOR_WEIGHT_ACCUMULATOR_HPP_

