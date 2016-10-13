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
#ifndef CLOTHO_BATCH_CROSSOVER_VECTOR_TASK_HPP_
#define CLOTHO_BATCH_CROSSOVER_VECTOR_TASK_HPP_

#include "clotho/data_spaces/population_space/population_space.hpp"
#include <boost/random/bernoulli_distribution.hpp>
#include "clotho/data_spaces/crossover/position_classifier.hpp"
#include "clotho/data_spaces/generators/position_distribution_helper.hpp"
#include "clotho/data_spaces/generators/crossover_event_distribution_helper.hpp"

#include "clotho/data_spaces/crossover/block_crossover.hpp"

namespace clotho {
namespace genetics {

template < class RNG, class MatePairType, class BlockType, class WeightType, class AlleleSpaceType >
class batch_crossover_task< RNG, MatePairType, population_space< BlockType, WeightType > , AlleleSpaceType > : public task {
public:
    typedef batch_crossover_task< RNG, MatePairType, population_space< BlockType, WeightType >, AlleleSpaceType > self_type;

    typedef population_space< BlockType, WeightType >   space_type;
    typedef AlleleSpaceType                             allele_type;
    typedef RNG                                         random_engine_type;

    typedef typename space_type::genome_type            genome_type;

    typedef typename space_type::individual_type        individual_type;

    typedef MatePairType                                        mate_pair_type;
    typedef typename mate_pair_type::iterator                   iterator;
    typedef typename mate_pair_type::const_iterator             const_iterator;

    typedef PositionClassifier< typename allele_type::position_vector >          classifier_type;
    typedef typename classifier_type::event_type                                event_type;

    typedef typename position_distribution_helper< typename allele_type::position_type >::type position_distribution_type;
    typedef typename crossover_event_distribution_helper< double >::type event_distribution_type;

    batch_crossover_task( random_engine_type * rng, space_type * parent, space_type * offspring, allele_type * alleles, unsigned int off_idx, const_iterator first, const_iterator last, double recomb_rate, double seq_bias ) :
        m_rng(rng)
        , m_parental( parent )
        , m_offspring( offspring )
        , m_alleles( alleles )
        , m_parents(first, last)
        , m_offspring_index(off_idx)
        , m_recomb_rate(recomb_rate)
        , m_seq_bias( seq_bias )
    { }

    batch_crossover_task( const self_type & other ) :
        m_rng( other.m_rng )
        , m_parental( other.m_parental )
        , m_offspring( other.m_offspring )
        , m_alleles( other.m_alleles )
        , m_parents( other.m_parents )
        , m_offspring_index( other.m_offspring_index )
        , m_recomb_rate( other.m_recomb_rate )
        , m_seq_bias( other.m_seq_bias )
    { }

    void operator()() {
        event_distribution_type     event_dist( m_recomb_rate);
        boost::random::bernoulli_distribution< double > bias_dist( m_seq_bias);

        iterator mate_it = m_parents.begin(), mate_end = m_parents.end();

#ifdef DEBUGGING
        BOOST_LOG_TRIVIAL(debug) << "Starting batch crossover task";
#endif // DEBUGGING
        size_t i = m_offspring_index;
        while( mate_it != mate_end ) {
            
            assert( mate_it->first < m_parental->individual_count() );

            unsigned int offs = mate_it->first;

            individual_type ind = m_parental->getIndividual( offs );
            genome_type c0 = m_offspring->create_sequence();

            event_type evts;

            fill_events( evts, event_dist( *m_rng ) );
            classifier_type cfier0( &m_alleles->getPositions(), evts );
            bool _swap = bias_dist( *m_rng );

#ifdef DEBUGGING
            BOOST_LOG_TRIVIAL(debug) << mate_it->first << "; Crossover s0: " << ind.first << " - " << ind.first.size() << "; s1: " << ind.second << " - " << ind.second.size() << "; child: " << c0.first << "; event size: " << evts.size();
#endif  // DEBUGGING
            run_crossover_task( cfier0, ind.first, ind.second, c0, _swap );

            assert( mate_it->second < m_parental->individual_count() );

            genome_type c1 = m_offspring->create_sequence();
            ind = m_parental->getIndividual( mate_it->second );
#ifdef DEBUGGING
            BOOST_LOG_TRIVIAL(debug) << mate_it->second << "; Crossover s0: " << ind.first << " - " << ind.first.size() << "; s1: " << ind.second << " - " << ind.second.size() << "; child: " << c1.first << "; event size: " << evts.size();
#endif  // DEBUGGING

            event_type evts1;

            fill_events( evts1, event_dist( *m_rng ) );

            classifier_type cfier1( &m_alleles->getPositions(), evts1 );

            _swap = bias_dist( *m_rng );
            run_crossover_task( cfier1, ind.first, ind.second, c1, _swap);

            m_offspring->setIndividual( i++, c0, c1 );

            ++mate_it;
        }

#ifdef DEBUGGING
        BOOST_LOG_TRIVIAL(debug) << "End Batch Crossover Task";
#endif // DEBUGGING
    }

    virtual ~batch_crossover_task() {}

protected:

    inline void fill_events( event_type & evt, unsigned int N ) {
        while( N-- ) {
            evt.push_back( m_pos_dist( *m_rng ) );
        }
    }

    event_type make_events( unsigned int N ) {
        event_type res;

        while( N-- ) {
            res.push_back( m_pos_dist( *m_rng ) );
        }

        return res;
    }

    void run_crossover_task( const classifier_type & cls, genome_type & top, genome_type & bottom, genome_type & res, bool should_swap_strands ) {

        if( should_swap_strands ) {
            run_crossover_task( cls, bottom, top, res );
        } else {
            run_crossover_task( cls, top, bottom, res );
        }
    }

    void run_crossover_task( const classifier_type & cls, const genome_type & top, const genome_type & bottom, genome_type & res ) {

        if( cls.event_count() == 0 ) {
            // there are no crossover cls
            // therefore, offspring strand will be a copy of the top strand
            res = top;
        } else {
            typedef typename space_type::base_genome_type::sequence_type::const_sequence_iterator const_iterator;
            typedef block_crossover< classifier_type, BlockType >       crossover_type;
            typedef BlockType                                           block_type;

            crossover_type xover( cls );
            
            const_iterator tb, te, bb, be;
            if( top ) {
                tb = top->begin_sequence();
                te = top->end_sequence();
            } else {
                te = tb;
            }

            if( bottom ) {
                bb = bottom->begin_sequence();
                be = bottom->end_sequence();
            } else {
                be = bb;
            }

            unsigned int i = 0;
            while( true ) {
                if( tb == te ) {
                    while( bb != be ) {
                        block_type o = xover.crossover( crossover_type::bit_helper_type::ALL_UNSET, *bb++, i++ );
                        res->append_sequence(o);
                    }
                    break;
                } else if( bb == be ) {
                    while( tb != te ) {
                        block_type o = xover.crossover( *tb++, crossover_type::bit_helper_type::ALL_UNSET, i++ );
                        res->append_sequence(o);
                    }
                    break;
                }

                block_type t = *tb++;
                block_type b = *bb++;

                block_type o = xover.crossover(t, b, i );

                res->append_sequence( o );
                ++i;
            }
        }
    }

    random_engine_type  * m_rng;
    space_type * m_parental, * m_offspring;
    allele_type         * m_alleles;

    mate_pair_type m_parents;

    unsigned int m_offspring_index;

    double m_recomb_rate, m_seq_bias;

    position_distribution_type  m_pos_dist;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_BATCH_CROSSOVER_VECTOR_TASK_HPP_

