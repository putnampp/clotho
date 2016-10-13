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
#ifndef ASSOCIATION_MATRIX_ALLELE_FREQUENCY_HPP_
#define ASSOCIATION_MATRIX_ALLELE_FREQUENCY_HPP_

#include "clotho/data_spaces/analysis/allele_frequency/allele_frequency_def.hpp"
#include "clotho/data_spaces/analysis/allele_frequency/frequency_evaluator.hpp"

#include "clotho/data_spaces/analysis/allele_frequency/margin_details.hpp"

#include "clotho/data_spaces/association_matrix.hpp"
//#include "clotho/utility/debruijn_bit_walker.hpp"

namespace clotho {
namespace genetics {

template < class BlockType, class AlignmentType, class SizeType >
class allele_frequency< association_matrix< BlockType, AlignmentType >, SizeType > : public margin_details< SizeType > {
public:
    typedef association_matrix< BlockType, AlignmentType >                  space_type;

    typedef typename margin_details< SizeType >::size_type                  size_type;
    typedef clotho::genetics::frequency_evaluator< space_type, size_type >        evaluator_type;

    allele_frequency() {}

    void evaluate( space_type & ss ) {
        this->resize( ss.column_count(), ss.row_count() );

        this->m_indices.reset();
        this->m_indices.flip();

        evaluator_type eval;

        eval( ss, this->m_indices, this->m_column_margin, this->m_row_margin );
    }

/**
 * Iterator sub population
 */
    template < class Iterator >
    void evaluate( space_type & ss, Iterator first, Iterator last ) {
        this->resize( ss.column_count(), ss.row_count() );

        this->m_indices.reset();

        size_type N = 0;
        while( first != last ) {
            size_type i = *first++;
            this->m_indices[i] = true;
            ++N;
        }

        evaluator_type eval;

        eval( ss, this->m_indices, this->m_column_margin, this->m_row_margin );
    }

    virtual ~allele_frequency() { }

protected:

//    void eval( space_type & ss ) {
//        m_column_margin.clear();
//
//        size_type M = ss.column_count();
//        size_type N = ss.block_column_count();
//        m_column_margin.resize( M );
//
//        size_type buffer[ bit_helper_type::BITS_PER_BLOCK ];
//        memset( buffer, 0, bit_helper_type::BITS_PER_BLOCK * sizeof( size_type ) );
//
//        typedef typename space_type::raw_block_pointer iterator;
//
//        size_type i = 0, j = M;
//        while( M ) {
//            iterator first = ss.begin_block_row( i );
//            iterator last = ss.end_block_row( i );
//
//            size_type k = 0;
//            while( first != last ) {
//                block_type b = *first++;
//                if( m_indices.test(k++) ) {
//                    while( b ) {
//                        unsigned int b_idx = bit_walker_type::unset_next_index( b );
//                        buffer[ b_idx ] += 1;
//                    }
//                }
//            }
//
//            size_type l = (( j < bit_helper_type::BITS_PER_BLOCK) ? j : bit_helper_type::BITS_PER_BLOCK);
//
//            memcpy( m_column_margin + i * bit_helper_type::BITS_PER_BLOCK, buffer, sizeof( size_type ) * l );
//            memset( buffer, 0, bit_helper_type::BITS_PER_BLOCK * sizeof( size_type ) );
//
//            j -= l;
//            ++i;
//        }
//    }

};

}   // namespace genetics
}   // namespace clotho

#endif  // ASSOCIATION_MATRIX_ALLELE_FREQUENCY_HPP_
