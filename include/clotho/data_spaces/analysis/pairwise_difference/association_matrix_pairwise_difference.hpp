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
#ifndef CLOTHO_ASSOCIATION_MATRIX_PAIRWISE_DIFFERENCE_HPP_
#define CLOTHO_ASSOCIATION_MATRIX_PAIRWISE_DIFFERENCE_HPP_

#include "clotho/data_spaces/analysis/pairwise_difference/pairwise_difference_def.hpp"
#include "clotho/data_spaces/analysis/pairwise_difference/pairwise_difference_evaluator.hpp"

#include "clotho/data_spaces/association_matrix.hpp"

#include <boost/dynamic_bitset.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/variance.hpp>

namespace clotho {
namespace genetics {

template < class BlockType, class AlignmentType >
class pairwise_difference< association_matrix< BlockType, AlignmentType > > {
public:

    typedef association_matrix< BlockType, AlignmentType >              space_type;
    typedef typename space_type::block_type                             block_type;

    typedef ac::accumulator_set< size_t, ac::stats< ac::tag::min, ac::tag::mean, ac::tag::max, ac::tag::variance, ac::tag::median > >                                                           result_type;
    typedef boost::dynamic_bitset<>                                     analyzed_set_type;

    typedef pairwise_difference_evaluator< BlockType, AlignmentType >   evaluator_type;
    pairwise_difference( ) { }

    void evaluate( space_type & ss ) {
        m_indices.resize( ss.row_count() );

        m_indices.reset();
        m_indices.flip();

        eval( ss );
    }

    template < class Iterator >
    void evaluate( space_type & ss, Iterator first, Iterator last ) {
        m_indices.resize( ss.row_count() );
        m_indices.reset();

        size_t N = 0;
        while( first != last ) {
            size_t i = *first++;
            m_indices[i] = true;
            ++N;
        }

        eval( ss );
    }
    
    result_type & getResults() {
        return m_results;
    }

    void recordResults( boost::property_tree::ptree & log ) {
        log.put( "min", ac::min( m_results ) );
        log.put( "max", ac::max( m_results ) );
        log.put( "mean", ac::mean( m_results ) );
        log.put( "median", ac::median( m_results ) );
        log.put( "variance", ac::variance( m_results ) );
    }

    virtual ~pairwise_difference() {}

protected:
    
    void eval( space_type & ss ) {
        typedef typename space_type::raw_block_pointer iterator;

        evaluator_type ee( ss.block_column_count() );

        size_t i = m_indices.find_first();
        while( i != analyzed_set_type::npos ) {
            size_t j = m_indices.find_next( i );
            while( j != analyzed_set_type::npos ) {

                iterator i_start = ss.begin_row( i );
                iterator i_end = ss.end_row(i);

                iterator j_start = ss.begin_row( j );
                iterator j_end = ss.end_row(j);

                size_t diff = ee( i_start, i_end, j_start, j_end );

                m_results( diff );

                j = m_indices.find_next( j );
            }

            i = m_indices.find_next(i);
        }
    }

    result_type         m_results;
    analyzed_set_type   m_indices;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_ASSOCIATION_MATRIX_PAIRWISE_DIFFERENCE_HPP_
