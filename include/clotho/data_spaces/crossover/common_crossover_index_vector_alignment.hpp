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
#ifndef CLOTHO_COMMON_CROSSOVER_INDEX_VECTOR_ALIGNMENT_HPP_
#define CLOTHO_COMMON_CROSSOVER_INDEX_VECTOR_ALIGNMENT_HPP_

#include "clotho/data_spaces/population_space/population_spaces.hpp"

#include <boost/random/bernoulli_distribution.hpp>
#include "clotho/data_spaces/crossover/position_classifier.hpp"
#include "clotho/data_spaces/generators/position_distribution_helper.hpp"
#include "clotho/data_spaces/generators/crossover_event_distribution_helper.hpp"

#include "clotho/data_spaces/crossover/index_crossover_method.hpp"

namespace clotho {
namespace genetics {

template < class RNG, class IndexType, class BlockType, class WeightType, class AlleleType >
class common_crossover< RNG, population_space< index_vector_alignment< IndexType, BlockType >, trait_space_vector< WeightType > >, AlleleType >  {
public:
    typedef RNG random_engine_type;

    typedef population_space< index_vector_alignment< IndexType, BlockType >, trait_space_vector< WeightType > >    space_type;
    typedef typename space_type::genome_pointer                                                             genome_pointer;
    typedef typename space_type::row_pointer                                                                row_pointer;

    typedef AlleleType                                                                                      allele_type;

    typedef PositionClassifier< typename allele_type::position_vector > classifier_type;
    typedef typename classifier_type::event_type                        event_type;

    typedef index_crossover_method< classifier_type, space_type >       method_type;

    typedef typename position_distribution_helper< typename allele_type::position_type >::type  position_distribution_type;
    typedef typename crossover_event_distribution_helper< double >::type                    event_distribution_type;

    common_crossover( random_engine_type * rng, allele_type * alleles, double recomb_rate, double bias_rate ) : 
        m_rng( rng )
        , m_alleles( alleles )
        , m_event_dist(recomb_rate)
        , m_bias_dist( bias_rate )
    { }

    void operator()( genome_pointer p0_start, genome_pointer p0_end, genome_pointer p1_start, genome_pointer p1_end, row_pointer offspring ) {
        event_type evts;
        fill_events( evts, m_event_dist(*m_rng));
        classifier_type cfier0( &m_alleles->getPositions(), evts );

        run_crossover_task( cfier0, p0_start, p0_end, p1_start, p1_end, offspring, m_bias_dist(*m_rng));
    }

    virtual ~common_crossover() {}

protected:

    void fill_events( event_type & evt, unsigned int N ) {
        while( N-- ) {
            evt.push_back( m_pos_dist( *m_rng ) );
        }
    }

    void run_crossover_task( const classifier_type & cls, genome_pointer p0_start, genome_pointer p0_end, genome_pointer p1_start, genome_pointer p1_end, row_pointer offspring, bool should_swap_strands ) {
        if( cls.event_count() == 0 ) {
            genome_pointer o_start = offspring->begin();
            size_t N = offspring->size();

            if( should_swap_strands ) {
                while( p1_start != p1_end ) {
//                    offspring->push_back( *p1_start++ );
                    if( N > 0 ) {
                        *o_start = *p1_start;
                        o_start++;
                        --N;
                    } else {
                        offspring->push_back( *p1_start );
                    }
                    ++p1_start;
                }
            } else {
                while( p0_start != p0_end ) {
//                    offspring->push_back( *p0_start++ );
                    if( N > 0 ) {
                        *o_start = *p0_start;
                        o_start++;
                        --N;
                    } else {
                        offspring->push_back( *p0_start );
                    }
                    p0_start++;
                }
            }

            while( N-- > 0 ) {
                offspring->pop_back();
            }
        } else if( should_swap_strands ) {
            method_type met(cls);
            met( p1_start, p1_end, p0_start, p0_end, offspring);
        } else {
            method_type met(cls);
            met( p0_start, p0_end, p1_start, p1_end, offspring);
        }   
    }

    random_engine_type          * m_rng;
    allele_type                 * m_alleles;
    event_distribution_type     m_event_dist;
    boost::random::bernoulli_distribution< double > m_bias_dist;
    position_distribution_type m_pos_dist;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_COMMON_CROSSOVER_INDEX_VECTOR_ALIGNMENT_HPP_
