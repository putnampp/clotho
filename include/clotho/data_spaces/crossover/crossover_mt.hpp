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
#ifndef CLOTHO_CROSSOVER_MT_GENERATOR_HPP_
#define CLOTHO_CROSSOVER_MT_GENERATOR_HPP_

//#include <boost/asio/io_service.hpp>
//#include <boost/bind.hpp>
//#include <boost/thread/thread.hpp>

#include <boost/random/bernoulli_distribution.hpp>

#include "clotho/data_spaces/crossover/crossover_task_list.hpp"
#include "clotho/recombination/sequence_bias_parameter.hpp"

#include "clotho/data_spaces/crossover/position_classifier.hpp"

#include "clotho/data_spaces/generators/position_distribution_helper.hpp"
#include "clotho/data_spaces/generators/crossover_event_distribution_helper.hpp"

#include "clotho/recombination/recombination_rate_parameter.hpp"

#include "clotho/data_spaces/task/thread_pool.hpp"

namespace clotho {
namespace genetics {

template < class RNG, class GeneticSpaceType >
class CrossoverMT {
public:
    typedef GeneticSpaceType    genetic_space_type;

    typedef typename genetic_space_type::allele_type        allele_type;
    typedef typename genetic_space_type::association_type   association_type;

    typedef typename association_type::sequence_vector      sequence_vector;
    typedef typename association_type::raw_block_pointer    sequence_iterator;
    typedef typename genetic_space_type::individual_id_type individual_id_type;

    typedef std::vector< std::pair< individual_id_type, individual_id_type > >  mate_pair_type;
    typedef typename mate_pair_type::iterator                                   iterator;

    typedef PositionClassifier< typename allele_type::position_type >          classifier_type;
    typedef typename classifier_type::event_type                                event_type;

    typedef typename position_distribution_helper< typename allele_type::position_type >::type position_distribution_type;
    typedef typename crossover_event_distribution_helper< double >::type event_distribution_type;

    CrossoverMT( RNG * rng, boost::property_tree::ptree & config ) :
        m_rng( rng )
    {
        sequence_bias_parameter< double > bias(config);
        m_bias_dist.param( typename boost::random::bernoulli_distribution< double >::param_type( bias.m_bias ) ); 

        recombination_rate_parameter< double > rho(config);
        m_event_dist.param( typename event_distribution_type::param_type( rho.m_rho ) );
    }

    template < class PoolType >
    void operator()( genetic_space_type * parental, mate_pair_type & parents, genetic_space_type * offspring, PoolType & pool ) {
        offspring->getSequenceSpace().clear();
        if( pool.pool_size() <= 1 ) {
            generate( parental, parents, offspring );
        } else {
            generate( parental, parents, offspring, pool );
        }
    }

    virtual ~CrossoverMT() {}

protected:
    
    void generate( genetic_space_type * parental, mate_pair_type & parents, genetic_space_type * offspring ) {
        iterator mate_it = parents.begin(), mate_end = parents.end();

        size_t i = 0;
        while( mate_it != mate_end ) {

            unsigned int idx = 2 * mate_it->first;
            sequence_vector s0 = parental->getSequenceSpace().getSequence( idx++ );
            sequence_vector s1 = parental->getSequenceSpace().getSequence( idx );

            sequence_iterator c = offspring->begin_sequence( i++ );

            event_type evts = make_events( );
            classifier_type cfier0( parental->getAlleleSpace().getPositions(), evts );
            bool _swap = m_bias_dist( *m_rng );
            make_crossover_tasks( cfier0, s0.first, s0.second, s1.first, s1.second, c, _swap ); 
            
            idx = 2 * mate_it->second;
            s0 = parental->getSequenceSpace().getSequence( idx++ );
            s1 = parental->getSequenceSpace().getSequence( idx );

            c = offspring->begin_sequence( i++ );

            event_type evts1 = make_events( );
            classifier_type cfier1( parental->getAlleleSpace().getPositions(), evts1 );
            _swap = m_bias_dist( *m_rng );
            make_crossover_tasks( cfier1, s0.first, s0.second, s1.first, s1.second, c, _swap);

            ++mate_it;
        }
    }

    template < class PoolType >
    void generate( genetic_space_type * parental, mate_pair_type & parents, genetic_space_type * offspring, PoolType & pool ) {
        
        iterator mate_it = parents.begin(), mate_end = parents.end();

        size_t i = 0;
        while( mate_it != mate_end ) {

            unsigned int idx = 2 * mate_it->first;
            sequence_vector s0 = parental->getSequenceSpace().getSequence( idx++ );
            sequence_vector s1 = parental->getSequenceSpace().getSequence( idx );

            sequence_iterator c = offspring->begin_sequence( i++ );

            event_type evts = make_events( );
            classifier_type cfier0( parental->getAlleleSpace().getPositions(), evts );
            bool _swap = m_bias_dist( *m_rng );
            make_crossover_tasks( cfier0, s0.first, s0.second, s1.first, s1.second, c, pool, _swap );
            
            idx = 2 * mate_it->second;
            s0 = parental->getSequenceSpace().getSequence( idx++ );
            s1 = parental->getSequenceSpace().getSequence( idx );

            c = offspring->begin_sequence( i++ );

            event_type evts1 = make_events( );
            classifier_type cfier1( parental->getAlleleSpace().getPositions(), evts1 );
            _swap = m_bias_dist( *m_rng );
            make_crossover_tasks( cfier1, s0.first, s0.second, s1.first, s1.second, c, pool, _swap );

            ++mate_it;
        }

        pool.sync();
    }

    event_type make_events( ) {
        event_type res;

        unsigned int N = m_event_dist( *m_rng );

        while( N-- ) {
            res.push_back( m_pos_dist( *m_rng ) );
        }

        return res;
    }

    RNG * m_rng;

    event_distribution_type     m_event_dist;
    position_distribution_type  m_pos_dist;
    boost::random::bernoulli_distribution< double > m_bias_dist;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_CROSSOVER_MT_GENERATOR_HPP_
