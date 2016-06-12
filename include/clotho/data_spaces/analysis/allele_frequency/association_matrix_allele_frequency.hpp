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

template < class BlockType, class AlignmentType >
class allele_frequency< association_matrix< BlockType, AlignmentType > > {
public:
    typedef association_matrix< BlockType, AlignmentType >             space_type;

    typedef clotho::genetics::frequency_evaluator< space_type >         evaluator_type;

    typedef size_t *                                                    result_type;
    typedef boost::dynamic_bitset<>                                     analyzed_set_type;

    allele_frequency() :
        m_results( NULL )
        , m_size(0)
        , m_allocated_size(0)
    {}

    void evaluate( space_type & ss ) {
        resize( ss );

        m_indices.reset();
        m_indices.flip();

        evaluator_type eval;

        eval( ss, m_indices, m_results );
    }

/**
 * Iterator sub population
 */
    template < class Iterator >
    void evaluate( space_type & ss, Iterator first, Iterator last ) {
        resize( ss );

        m_indices.reset();

        size_t N = 0;
        while( first != last ) {
            size_t i = *first++;
            m_indices[i] = true;
            ++N;
        }

        evaluator_type eval;

        eval( ss, m_indices, m_results );
    }

    std::shared_ptr< std::vector< size_t > > getResults() {
        if( m_results == NULL || m_size == 0 ) {
            return std::shared_ptr< std::vector< size_t > >();
        }
        std::shared_ptr< std::vector< size_t > > t( new std::vector< size_t >( m_results, m_results + m_size) );
        return t;
    }

    void recordResults( boost::property_tree::ptree & log ) {
                    
        result_type b = m_results, e = m_results + m_size;
 
        typedef ac::accumulator_set< size_t, ac::stats< ac::tag::min, ac::tag::mean, ac::tag::max, ac::tag::variance, ac::tag::median, ac::tag::count > > accumulator_t;

        accumulator_t   accum;       

        boost::property_tree::ptree dist, lfreq;

        std::vector< size_t > freq( m_indices.size(), 0 );
        
        while( b != e ) {
            size_t v = *b++;
            clotho::utility::add_value_array( dist, v );
            accum( v );
            freq[v]++;
        }

        clotho::utility::add_value_array(lfreq, freq.begin(), freq.end() );

        log.put_child( "distribution", dist );
        log.put_child( "frequency_distribution", lfreq );

        log.put( "min", ac::min( accum ) );
        log.put( "max", ac::max( accum ) );
        log.put( "mean", ac::mean( accum ) );
        log.put( "median", ac::median( accum ) );
        log.put( "variance", ac::variance( accum ) );
        log.put( "total", ac::count(accum) );
    }

    virtual ~allele_frequency() {
        if( m_results != NULL ) {
            delete [] m_results;
        }
    }

protected:

//    void eval( space_type & ss ) {
//        m_results.clear();
//
//        size_t M = ss.column_count();
//        size_t N = ss.block_column_count();
//        m_results.resize( M );
//
//        size_t buffer[ bit_helper_type::BITS_PER_BLOCK ];
//        memset( buffer, 0, bit_helper_type::BITS_PER_BLOCK * sizeof( size_t ) );
//
//        typedef typename space_type::raw_block_pointer iterator;
//
//        size_t i = 0, j = M;
//        while( M ) {
//            iterator first = ss.begin_block_row( i );
//            iterator last = ss.end_block_row( i );
//
//            size_t k = 0;
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
//            size_t l = (( j < bit_helper_type::BITS_PER_BLOCK) ? j : bit_helper_type::BITS_PER_BLOCK);
//
//            memcpy( m_results + i * bit_helper_type::BITS_PER_BLOCK, buffer, sizeof( size_t ) * l );
//            memset( buffer, 0, bit_helper_type::BITS_PER_BLOCK * sizeof( size_t ) );
//
//            j -= l;
//            ++i;
//        }
//    }

    void resize( space_type & ss ) {
        size_t M = ss.column_count();

        if( M > m_allocated_size ) {
            if( m_results != NULL ) {
                delete [] m_results;
            }

            m_results = new size_t[ M ];

            m_allocated_size = M;
        }

        m_size = M;

        memset( m_results, 0, sizeof(size_t) * m_allocated_size );

        m_indices.resize( ss.row_count() );

    }

    result_type         m_results;
    analyzed_set_type   m_indices;

    size_t              m_size, m_allocated_size;
};

}   // namespace genetics
}   // namespace clotho

#endif  // ASSOCIATION_MATRIX_ALLELE_FREQUENCY_HPP_
