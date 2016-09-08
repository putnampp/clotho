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

#include "clotho/data_spaces/association_matrix.hpp"
#include "clotho/utility/debruijn_bit_walker.hpp"

#include <boost/dynamic_bitset.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/count.hpp>

namespace ac = boost::accumulators;

namespace clotho {
namespace genetics {

template < class BlockType, class AlignmentType, class SizeType >
class allele_frequency< association_matrix< BlockType, AlignmentType >, SizeType > {
public:
    typedef association_matrix< BlockType, AlignmentType >             space_type;

    typedef SizeType                                                      size_type;
    typedef size_type *                                                   result_type;
    typedef clotho::genetics::frequency_evaluator< space_type, size_type >        evaluator_type;

    typedef boost::dynamic_bitset<>                                     row_subset_type;

    allele_frequency() :
        m_column_margin( NULL )
        , m_row_margin( NULL )
        , m_column_margin_size(0)
        , m_column_allocated_size(0)
        , m_row_margin_size(0)
        , m_row_allocated_size(0)
    {}

    void evaluate( space_type & ss ) {
        resize( ss );

        m_indices.reset();
        m_indices.flip();

        evaluator_type eval;

        eval( ss, m_indices, m_column_margin, m_row_margin );
    }

/**
 * Iterator sub population
 */
    template < class Iterator >
    void evaluate( space_type & ss, Iterator first, Iterator last ) {
        resize( ss );

        m_indices.reset();

        size_type N = 0;
        while( first != last ) {
            size_type i = *first++;
            m_indices[i] = true;
            ++N;
        }

        evaluator_type eval;

        eval( ss, m_indices, m_column_margin, m_row_margin );
    }

    std::shared_ptr< std::vector< size_type > > getResults() {
        if( m_column_margin == NULL || m_column_margin_size == 0 ) {
            return std::shared_ptr< std::vector< size_type > >();
        }
        std::shared_ptr< std::vector< size_type > > t( new std::vector< size_type >( m_column_margin, m_column_margin + m_column_margin_size) );
        return t;
    }

    void recordResults( boost::property_tree::ptree & log ) {
 
        typedef ac::accumulator_set< size_type, ac::stats< ac::tag::min, ac::tag::mean, ac::tag::max, ac::tag::variance, ac::tag::median, ac::tag::count > > accumulator_t;

        accumulator_t   col_accum, row_accum;

        boost::property_tree::ptree dist, lfreq;

        // maximum frequency for an allele is bounded by the number of sequences (rows) in the population
        std::vector< size_type > freq( m_indices.size(), 0 );

        // evaluate allele (column) statistics
        size_type i = 0;
        while( i < m_column_margin_size ) {
            size_type v = m_column_margin[i++];
            clotho::utility::add_value_array( dist, v );
            col_accum( v );
            freq[v]++;
        }


        clotho::utility::add_value_array(lfreq, freq.begin(), freq.end() );

        log.put_child( "distribution", dist );
        log.put_child( "frequency_distribution", lfreq );

        log.put( "stats.sequences_per_allele.min", ac::min( col_accum ) );
        log.put( "stats.sequences_per_allele.max", ac::max( col_accum ) );
        log.put( "stats.sequences_per_allele.mean", ac::mean( col_accum ) );
        log.put( "stats.sequences_per_allele.median", ac::median( col_accum ) );
        log.put( "stats.sequences_per_allele.variance", ac::variance( col_accum ) );
        log.put( "stats.sequences_per_allele.total", ac::count(col_accum) );

        // evaluate sequence (row) statistics
        i = 0;
        while( i < m_row_margin_size ) {
            if( m_indices.test(i) ) {
                row_accum( m_row_margin[ i ] );
            }
            ++i;
        }

        log.put( "stats.alleles_per_sequence.min", ac::min( row_accum ) );
        log.put( "stats.alleles_per_sequence.max", ac::max( row_accum ) );
        log.put( "stats.alleles_per_sequence.mean", ac::mean( row_accum ) );
        log.put( "stats.alleles_per_sequence.median", ac::median( row_accum ) );
        log.put( "stats.alleles_per_sequence.variance", ac::variance( row_accum ) );
        log.put( "stats.alleles_per_sequence.total", ac::count(row_accum) );
    }

    virtual ~allele_frequency() {
        if( m_column_margin != NULL ) {
            delete [] m_column_margin;
        }
    }

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

/**
 *
 * Grows vectors if necessary
 * Always clears the count vectors to be all zeros
 *
 */
    void resize( space_type & ss ) {
        size_type M = ss.column_count();

        if( M > m_column_allocated_size ) {
            if( m_column_margin != NULL ) {
                delete [] m_column_margin;
            }

            m_column_margin = new size_type[ M ];

            m_column_allocated_size = M;
        }

        m_column_margin_size = M;

        memset( m_column_margin, 0, sizeof(size_type) * m_column_allocated_size );

        m_indices.resize( ss.row_count() );
        if( ss.row_count() > m_row_allocated_size ) {
            if( m_row_margin != NULL ) {
                delete [] m_row_margin;
            }

            m_row_margin = new size_type[ ss.row_count() ];
            m_row_allocated_size = ss.row_count();
        }

        m_row_margin_size = ss.row_count();

        memset( m_row_margin, 0, sizeof(size_type) * m_row_allocated_size );
    }

    result_type         m_column_margin;
    result_type         m_row_margin;
    row_subset_type     m_indices;

    size_type              m_column_margin_size, m_column_allocated_size;
    size_type              m_row_margin_size, m_row_allocated_size;
};

}   // namespace genetics
}   // namespace clotho

#endif  // ASSOCIATION_MATRIX_ALLELE_FREQUENCY_HPP_
