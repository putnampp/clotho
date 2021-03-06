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
#ifndef CLOTHO_BATCH_CROSSOVER_TASKWE_HPP_
#define CLOTHO_BATCH_CROSSOVER_TASKWE_HPP_

#include "clotho/data_spaces/task/task.hpp"
#include "clotho/data_spaces/task/thread_pool.hpp"

#include <boost/random/bernoulli_distribution.hpp>
#include "clotho/data_spaces/crossover/position_classifier.hpp"
#include "clotho/data_spaces/generators/position_distribution_helper.hpp"
#include "clotho/data_spaces/generators/crossover_event_distribution_helper.hpp"
#include "clotho/data_spaces/crossover/crossover_task_list.hpp"

namespace clotho {
namespace genetics {

template <  class MatePairType, class SequenceSpaceType, class AlleleSpaceType >
class batch_crossover_taskwe;

/* : public task {
public:
    typedef batch_crossover_taskwe< MatePairType, SequenceSpaceType, AlleleSpaceType >   self_type;

    typedef SequenceSpaceType                               sequence_space_type;
    typedef AlleleSpaceType                                 allele_type;
    typedef RNG                                             random_engine_type;

    typedef typename sequence_space_type::sequence_vector      sequence_vector;
    typedef typename sequence_space_type::raw_block_pointer    sequence_iterator;
//    typedef typename sequence_space_type::individual_id_type   individual_id_type;

//    typedef std::vector< std::pair< individual_id_type, individual_id_type > >  mate_pair_type;
    typedef MatePairType                                        mate_pair_type;
    typedef typename mate_pair_type::iterator                   iterator;
    typedef typename mate_pair_type::const_iterator             const_iterator;

    typedef PositionClassifier< typename allele_type::position_vector >          classifier_type;
    typedef typename classifier_type::event_type                                event_type;

    typedef typename position_distribution_helper< typename allele_type::position_type >::type position_distribution_type;
    typedef typename crossover_event_distribution_helper< double >::type event_distribution_type;

    batch_crossover_taskwe ( random_engine_type * rng, sequence_space_type * parent, sequence_space_type * offspring, allele_type * alleles, unsigned int off_idx, const_iterator first, const_iterator last, double recomb_rate, double seq_bias) :
        m_rng(rng)
        , m_parental( parent )
        , m_offspring( offspring )
        , m_alleles( alleles )
        , m_parents(first, last )
        , m_offspring_index(off_idx)
        , m_recomb_rate(recomb_rate)
        , m_seq_bias( seq_bias )
    { }

    batch_crossover_taskwe( const self_type & other ) :
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
        size_t i = 2 * m_offspring_index;
        while( mate_it != mate_end ) {

            unsigned int idx = 2 * mate_it->first;

            assert( idx < m_parental->row_count() );

            sequence_vector s0 = m_parental->getSequence( idx++ );
            sequence_vector s1 = m_parental->getSequence( idx );

            sequence_vector c = m_offspring->getSequence( i++ );

            event_type evts, evts1;

            fill_events( evts, event_dist( *m_rng ) );
            classifier_type cfier0( &m_alleles->getPositions(), evts );
            bool _swap = bias_dist( *m_rng );

#ifdef DEBUGGING
            BOOST_LOG_TRIVIAL(debug) << mate_it->first << "; Crossover s0: " << s0.first << " - " << s0.second << "; s1: " << s1.first << " - " << s1.second << "; child: " << c.first << "; event size: " << evts.size();
#endif  // DEBUGGING
            run_crossover_task( cfier0, s0.first, s0.second, s1.first, s1.second, c.first, _swap ); 
            
            idx = 2 * mate_it->second;
            assert( idx < m_parental->row_count() );

            s0 = m_parental->getSequence( idx++ );
            s1 = m_parental->getSequence( idx );

            c = m_offspring->getSequence( i++ );
#ifdef DEBUGGING
            BOOST_LOG_TRIVIAL(debug) << mate_it->second << "; Crossover s0: " << s0.first << " - " << s0.second << "; s1: " << s1.first << " - " << s1.second << "; child: " << c.first << "; event size: " << evts.size();
#endif  // DEBUGGING

            fill_events( evts1, event_dist( *m_rng ) );

            classifier_type cfier1( &m_alleles->getPositions(), evts1 );
            _swap = bias_dist( *m_rng );
            run_crossover_task( cfier1, s0.first, s0.second, s1.first, s1.second, c.first, _swap);

            ++mate_it;
        }

#ifdef DEBUGGING
        BOOST_LOG_TRIVIAL(debug) << "End Batch Crossover Task";
#endif // DEBUGGING
    }

    virtual ~batch_crossover_taskwe() {}

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

    template < class Classifier, class BlockType >
    void run_crossover_task( const Classifier & cls, BlockType * s0, unsigned int s0_len, BlockType * s1, unsigned int s1_len, BlockType * res ) {
        typedef copy_crossover_task< BlockType > copy_type;
        typedef segment_crossover_task< Classifier, BlockType > segment_type;
        typedef tail_crossover_task< Classifier, BlockType, top_strand_tail > top_tail_type;
        typedef tail_crossover_task< Classifier, BlockType, bottom_strand_tail > bottom_tail_type;

        if( s0_len == 0 && s1_len == 0 ) return;

        if( cls.event_count() == 0 ) {
            // there are no crossover cls
            // therefore, offspring strand will be a copy of the top strand
            copy_type t( s0, res, s0_len );
            t();
        } else if( s0_len < s1_len ) {
            // top strand is shorter than bottom strand
            segment_type s( cls, s0, s1, res, 0, s0_len );
            s();

            bottom_tail_type t( cls, s1, res, s0_len, s1_len );
            t();
        } else if( s1_len < s0_len ) {
            // bottom strand is shorter than top strand
            segment_type s( cls, s0, s1, res, 0, s1_len );
            s();

            top_tail_type t( cls, s0, res, s1_len, s0_len );
            t();
        } else {
            // both strands are equal in length
            segment_type s( cls, s0, s1, res, 0, s0_len );
            s();
        }
    }

    template < class Classifier, class BlockType >
    void run_crossover_task( const Classifier & cls, BlockType * s0, unsigned int s0_len, BlockType * s1, unsigned int s1_len, BlockType * res, bool should_swap_strands ) {

        if( should_swap_strands ) {
            run_crossover_task( cls, s1, s1_len, s0, s0_len, res );
        } else {
            run_crossover_task( cls, s0, s0_len, s1, s1_len, res );
        }
    }

    random_engine_type  * m_rng;
    sequence_space_type * m_parental, * m_offspring;
    allele_type         * m_alleles;

    mate_pair_type m_parents;

    unsigned int m_offspring_index;

    double m_recomb_rate, m_seq_bias;

    position_distribution_type  m_pos_dist;
};
*/

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_BATCH_CROSSOVER_TASKWE_HPP_

