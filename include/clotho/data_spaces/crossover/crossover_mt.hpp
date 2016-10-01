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

template < class RNG, class SequenceSpaceType, class AlleleSpaceType >
class CrossoverMT {
public:
    typedef SequenceSpaceType   sequence_space_type;
    typedef AlleleSpaceType     allele_space_type;

    typedef typename sequence_space_type::sequence_vector      sequence_vector;
    typedef typename sequence_space_type::raw_block_pointer    sequence_iterator;
//    typedef typename sequence_space_type::individual_id_type individual_id_type;
//
//    typedef std::vector< std::pair< individual_id_type, individual_id_type > >  mate_pair_type;
//    typedef typename mate_pair_type::iterator                                   iterator;

    typedef PositionClassifier< typename allele_space_type::position_vector >   classifier_type;
    typedef typename classifier_type::event_type                                event_type;

    typedef typename position_distribution_helper< typename allele_space_type::position_type >::type position_distribution_type;
    typedef typename crossover_event_distribution_helper< double >::type event_distribution_type;

    CrossoverMT( RNG * rng, boost::property_tree::ptree & config ) :
        m_rng( rng )
    {
        sequence_bias_parameter< double > bias(config);
        m_bias_dist.param( typename boost::random::bernoulli_distribution< double >::param_type( bias.m_bias ) ); 

        recombination_rate_parameter< double > rho(config);
        m_event_dist.param( typename event_distribution_type::param_type( rho.m_rho ) );
    }

    template < class SelectionType, class PoolType >
    void operator()( const SelectionType & parents, sequence_space_type * parental, sequence_space_type * offspring, allele_space_type * alleles, PoolType & pool ) {
        offspring->clear();
        if( pool.pool_size() <= 1 ) {
            generate( parents, parental, offspring, alleles );
        } else {
            generate( parents, parental, offspring, alleles, pool );
        }
    }

    virtual ~CrossoverMT() {}

protected:
    
    template < class SelectionType >
    void generate( const SelectionType & parents, sequence_space_type * parental, sequence_space_type * offspring, allele_space_type * alleles ) {
        typedef typename SelectionType::const_iterator iterator;
        iterator mate_it = parents.begin(), mate_end = parents.end();

        size_t i = 0;
        while( mate_it != mate_end ) {

            unsigned int idx = 2 * mate_it->first;
            sequence_vector s0 = parental->getSequence( idx++ );
            sequence_vector s1 = parental->getSequence( idx );

            sequence_vector c = offspring->getSequence( i++ );

            event_type evts = make_events( );
            classifier_type cfier0( &alleles->getPositions(), evts );
            bool _swap = m_bias_dist( *m_rng );
            make_crossover_tasks( cfier0, s0.first, s0.second, s1.first, s1.second, c.first, _swap ); 
            
            idx = 2 * mate_it->second;
            s0 = parental->getSequence( idx++ );
            s1 = parental->getSequence( idx );

            c = offspring->getSequence( i++ );

            event_type evts1 = make_events( );
            classifier_type cfier1( &alleles->getPositions(), evts1 );
            _swap = m_bias_dist( *m_rng );
            make_crossover_tasks( cfier1, s0.first, s0.second, s1.first, s1.second, c.first, _swap);

            ++mate_it;
        }
    }

    template < class SelectionType, class PoolType >
    void generate( const SelectionType & parents, sequence_space_type * parental, sequence_space_type * offspring, allele_space_type * alleles, PoolType & pool ) {
        typedef typename SelectionType::const_iterator iterator;
        
        iterator mate_it = parents.begin(), mate_end = parents.end();

        size_t i = 0;
        while( mate_it != mate_end ) {

            unsigned int idx = 2 * mate_it->first;
            sequence_vector s0 = parental->getSequence( idx++ );
            sequence_vector s1 = parental->getSequence( idx );

            sequence_vector c = offspring->getSequence( i++ );

            event_type evts = make_events( );
            classifier_type cfier0( &alleles->getPositions(), evts );
            bool _swap = m_bias_dist( *m_rng );
            make_crossover_tasks( cfier0, s0.first, s0.second, s1.first, s1.second, c.first, pool, _swap );
            
            idx = 2 * mate_it->second;
            s0 = parental->getSequence( idx++ );
            s1 = parental->getSequence( idx );

            c = offspring->getSequence( i++ );

            event_type evts1 = make_events( );
            classifier_type cfier1( &alleles->getPositions(), evts1 );
            _swap = m_bias_dist( *m_rng );
            make_crossover_tasks( cfier1, s0.first, s0.second, s1.first, s1.second, c.first, pool, _swap );

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
