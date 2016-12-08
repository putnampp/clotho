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
#ifndef CLOTHO_BATCH_CROSSOVER_MTWE_GENERATOR_HPP_
#define CLOTHO_BATCH_CROSSOVER_MTWE_GENERATOR_HPP_

#include "clotho/data_spaces/crossover/batch_crossover_tasks.hpp"
#include "clotho/recombination/sequence_bias_parameter.hpp"
#include "clotho/recombination/recombination_rate_parameter.hpp"
#include "clotho/data_spaces/task/thread_pool.hpp"

#include "clotho/data_spaces/generators/position_distribution_helper.hpp"
#include "clotho/data_spaces/generators/crossover_event_distribution_helper.hpp"

#include <boost/random/bernoulli_distribution.hpp>

namespace clotho {
namespace genetics {

template < class RNG, class SequenceSpaceType, class AlleleSpaceType >
class BatchCrossoverMTWE {
public:
    typedef SequenceSpaceType       sequence_space_type;
    typedef AlleleSpaceType         allele_space_type;

    typedef typename position_distribution_helper< typename allele_space_type::position_type >::type  position_distribution_type;
    typedef typename crossover_event_distribution_helper< double >::type                    event_distribution_type;

    BatchCrossoverMTWE( RNG * rng, boost::property_tree::ptree & config ) :
        m_rng( rng )
        , m_seq_bias(config)
        , m_recomb_rate(config)
    { }

    template < class SelectionType, class PoolType >
    void operator()( const SelectionType & parents, sequence_space_type * parental, sequence_space_type * offspring, allele_space_type * alleles, PoolType & pool ) {

        batch_generate( parents, parental, offspring, alleles, pool );
    }

    virtual ~BatchCrossoverMTWE() {}

protected:
    
    template < class SelectionType, class PoolType >
    void batch_generate( const SelectionType & parents, sequence_space_type * parental, sequence_space_type * offspring, allele_space_type * alleles, PoolType & pool ) {
#ifdef  DEBUGGING
        BOOST_LOG_TRIVIAL(debug) << "Launching crossover batch jobs";
#endif  // DEBUGGING
        typedef typename SelectionType::mate_pair_vector                            mate_pair_type;
        typedef batch_crossover_taskwe< mate_pair_type, SequenceSpaceType, AlleleSpaceType >     crossover_task_type;

        typedef typename crossover_task_type::event_pool_type                       event_pool_type;
        typedef typename crossover_task_type::event_lookup_type                     event_lookup_type;
        typedef typename crossover_task_type::bias_pool_type                        bias_pool_type;

        event_pool_type ep = crossover_task_type::make_event_pool();
        event_lookup_type el;

        fill_event_pool( ep, el, 2 * parents.size() );

        bias_pool_type bp;
        fill_bias_pool( bp, 2 * parents.size() );
    
        const size_t TC = pool.pool_size() + 1;   // + 1 for master thread

        const size_t  BATCH_SIZE = (parents.size() / TC) + ((parents.size() % TC > 0)? 1 : 0);


        unsigned int off_idx = 0;
        while( off_idx + BATCH_SIZE < parents.size() ) {
            unsigned int off_end = off_idx + BATCH_SIZE;

#ifdef DEBUGGING
            BOOST_LOG_TRIVIAL(info) << "Batch Crossover: [" << off_idx << ", " << off_end << ");";
#endif  // DEBUGGING

            crossover_task_type x( parental, offspring, alleles, ep, off_idx, parents.begin() + off_idx, parents.begin() + off_end, el.begin() + 2 * off_idx, el.begin() + 2 * off_end + 1,  bp.begin() + 2 * off_idx, bp.begin() + 2 * off_end );
            pool.post( x );

            off_idx = off_end;
        }

        if( off_idx < parents.size() ){
#ifdef DEBUGGING
            BOOST_LOG_TRIVIAL(info) << "Batch Crossover: [" << off_idx << ", " << parents.size() << ");";
#endif  // DEBUGGING

            crossover_task_type t( parental, offspring, alleles, ep, off_idx, parents.begin() + off_idx, parents.end(), el.begin() + 2 * off_idx, el.end(), bp.begin() + 2 * off_idx, bp.end());
            t();
        }

        pool.sync();
#ifdef DEBUGGING
        BOOST_LOG_TRIVIAL(debug) << "Thread pool synced";
#endif // DEBUGGING
    }

    template < class EventPool, class EventLookup >
    void fill_event_pool( EventPool & events, EventLookup & lookup, unsigned int offspring_pop_size ) {
        event_distribution_type event_dist( m_recomb_rate.m_rho );
        position_distribution_type pos_dist;

        while( offspring_pop_size-- ) {
            lookup.push_back( events->size() );
            unsigned int N = event_dist( *m_rng );
            while( N-- ) {
                events->push_back( pos_dist( *m_rng ) );
            }
        }

        lookup.push_back( events->size() );
    }

    template < class BiasPool >
    void fill_bias_pool( BiasPool & bp, unsigned int offspring_pop_size ) {
        boost::random::bernoulli_distribution< double > bias_dist( m_seq_bias.m_bias );

        while( offspring_pop_size-- ) {
            bp.push_back( bias_dist( *m_rng ) );
        }
    }

    RNG * m_rng;
    sequence_bias_parameter< double > m_seq_bias;
    recombination_rate_parameter< double > m_recomb_rate;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_BATCH_CROSSOVER_MTWE_GENERATOR_HPP_
