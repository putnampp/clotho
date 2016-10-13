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
#ifndef CLOTHO_MARGIN_DETAILS_HPP_
#define CLOTHO_MARGIN_DETAILS_HPP_

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

template < class SizeType >
class margin_details {
public:

    typedef SizeType size_type;
    typedef size_type * result_type;

    typedef boost::dynamic_bitset<>                                     row_subset_type;

    margin_details() :
        m_column_margin( NULL )
        , m_row_margin( NULL )
        , m_column_margin_size(0)
        , m_column_allocated_size(0)
        , m_row_margin_size(0)
        , m_row_allocated_size(0)
    {}

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

    std::shared_ptr< std::vector< size_type > > getResults() {
        if( m_column_margin == NULL || m_column_margin_size == 0 ) {
            return std::shared_ptr< std::vector< size_type > >();
        }
        std::shared_ptr< std::vector< size_type > > t( new std::vector< size_type >( m_column_margin, m_column_margin + m_column_margin_size) );
        return t;
    }

    virtual ~margin_details() {
        if( m_column_margin != NULL ) {
            delete [] m_column_margin;
        }

        if( m_row_margin != NULL ) {
            delete [] m_row_margin;
        }
    }

protected:
/**
 *
 * Grows vectors if necessary
 * Always clears the count vectors to be all zeros
 *
 */
    void resize( size_type M, size_type N ) {

        if( M > m_column_allocated_size ) {
            if( m_column_margin != NULL ) {
                delete [] m_column_margin;
            }

            m_column_margin = new size_type[ M ];

            m_column_allocated_size = M;
        }

        m_column_margin_size = M;

        memset( m_column_margin, 0, sizeof(size_type) * m_column_allocated_size );

        m_indices.resize( N );
        if( N > m_row_allocated_size ) {
            if( m_row_margin != NULL ) {
                delete [] m_row_margin;
            }

            m_row_margin = new size_type[ N ];
            m_row_allocated_size = N;
        }

        m_row_margin_size = N;

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
#endif  // CLOTHO_MARGIN_DETAILS_HPP_
