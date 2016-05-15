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

namespace ac = boost::accumulators;

namespace clotho {
namespace genetics {

template < class BlockType >
class allele_frequency< association_matrix< BlockType, column_aligned > > {
public:
    typedef association_matrix< BlockType, column_aligned >             space_type;

    typedef typename space_type::block_type                             block_type;
    typedef typename space_type::bit_helper_type                        bit_helper_type;
    typedef typename clotho::utility::debruijn_bit_walker< block_type > bit_walker_type;

    typedef std::vector< size_t >                                       result_type;
    typedef boost::dynamic_bitset<>                                     analyzed_set_type;

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

    result_type &   getResults() {
        return m_results;
    }

    void recordResults( boost::property_tree::ptree & log ) {
        typedef typename result_type::iterator          iterator;
                    
        iterator b = m_results.begin(), e = m_results.end();
 
        typedef ac::accumulator_set< size_t, ac::stats< ac::tag::min, ac::tag::mean, ac::tag::max, ac::tag::variance, ac::tag::median > > accumulator_t;

        accumulator_t   accum;       

        boost::property_tree::ptree dist;
        
        while( b != e ) {
            size_t v = *b++;
            clotho::utility::add_value_array( dist, v );
            accum( v );
        }

        log.put_child( "distribution", dist );

        log.put( "min", ac::min( accum ) );
        log.put( "max", ac::max( accum ) );
        log.put( "mean", ac::mean( accum ) );
        log.put( "median", ac::median( accum ) );
        log.put( "variance", ac::variance( accum ) );
    }

    virtual ~allele_frequency() {}

protected:

    void eval( space_type & ss ) {
        m_results.clear();

        size_t M = ss.column_count();
        size_t N = ss.block_column_count();
        m_results.resize( M );

        size_t buffer[ bit_helper_type::BITS_PER_BLOCK ];
        memset( buffer, 0, bit_helper_type::BITS_PER_BLOCK * sizeof( size_t ) );

        typedef typename space_type::block_iterator    iterator;
        iterator b_it = ss.raw_iterator();

        size_t i = 0, j = 0;
        while( b_it.hasNext() ) {

            block_type b = b_it.next();

            if( m_indices[ i ] ) {
                while( b ) {
                    unsigned int b_idx = bit_walker_type::unset_next_index( b );
                    buffer[ b_idx ] += 1;
                }
            }

            if( ++i >= N ) {
                size_t k = 0;
                size_t l = (((j + bit_helper_type::BITS_PER_BLOCK ) > M) ? (M - j) : bit_helper_type::BITS_PER_BLOCK);
                size_t m = (l & (size_t)3); // l % 4

                while( k < m ) {
                    m_results[ j + k ] = buffer[ k ];
                    ++k;
                }
                
                while( k < l ) {
                    m_results[ j + k ] = buffer[ k ];
                    m_results[ j + k + 1 ] = buffer[ k + 1 ];
                    m_results[ j + k + 2 ] = buffer[ k + 2 ];
                    m_results[ j + k + 3 ] = buffer[ k + 3 ];
                    k += 4;
                }

                memset( buffer, 0, bit_helper_type::BITS_PER_BLOCK * sizeof( size_t ) );

                i = 0;
                j += bit_helper_type::BITS_PER_BLOCK;
            }
        }
    }

    result_type         m_results;
    analyzed_set_type   m_indices;
};

}   // namespace genetics
}   // namespace clotho

#endif  // ASSOCIATION_MATRIX_ALLELE_FREQUENCY_HPP_
