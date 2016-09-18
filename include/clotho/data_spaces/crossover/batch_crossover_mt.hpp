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
#ifndef CLOTHO_BATCH_CROSSOVER_MT_GENERATOR_HPP_
#define CLOTHO_BATCH_CROSSOVER_MT_GENERATOR_HPP_

#include "clotho/data_spaces/crossover/batch_crossover_task.hpp"
#include "clotho/recombination/sequence_bias_parameter.hpp"
#include "clotho/recombination/recombination_rate_parameter.hpp"
#include "clotho/data_spaces/task/thread_pool.hpp"

#include "clotho/data_spaces/crossover/crossover_task_list.hpp"

namespace clotho {
namespace genetics {

template < class RNG, class GeneticSpaceType >
class BatchCrossoverMT {
public:
    typedef GeneticSpaceType    genetic_space_type;

    typedef typename genetic_space_type::individual_id_type individual_id_type;

    typedef std::vector< std::pair< individual_id_type, individual_id_type > >  mate_pair_type;
    typedef batch_crossover_task< RNG, GeneticSpaceType >                       crossover_task_type;

    typedef clotho::genetics::thread_pool< RNG >                                pool_type;

    BatchCrossoverMT( RNG * rng, boost::property_tree::ptree & config ) :
        m_rng( rng )
        , m_seq_bias(config)
        , m_recomb_rate(config)
    { }

    void operator()( genetic_space_type * parental, mate_pair_type & parents, genetic_space_type * offspring, pool_type & pool ) {
        offspring->getSequenceSpace().clear();

        batch_generate( parental, parents, offspring, pool );
    }

    virtual ~BatchCrossoverMT() {}

protected:
    
    void batch_generate( genetic_space_type * parental, mate_pair_type & parents, genetic_space_type * offspring, pool_type & pool ) {
#ifdef  DEBUGGING
        BOOST_LOG_TRIVIAL(debug) << "Launching crossover batch jobs";
#endif  // DEBUGGING
        size_t tc = pool.pool_size() + 1;   // + 1 for master thread

        size_t  batch_size = parents.size() / tc;
        batch_size += ((parents.size() % tc > 0)? 1 : 0);

        unsigned int off_idx = 0, j = 0;
        while( off_idx + batch_size < parents.size() ) {
            unsigned int off_end = off_idx + batch_size;

            pool.post( crossover_task_type( pool.getRNG(++j), parental, offspring, off_idx, parents.begin() + off_idx, parents.begin() + off_end, m_recomb_rate.m_rho, m_seq_bias.m_bias ));

            off_idx = off_end;
        }

        crossover_task_type t( m_rng, parental, offspring, off_idx, parents.begin() + off_idx, parents.end(), m_recomb_rate.m_rho, m_seq_bias.m_bias);
        t();

        pool.sync();
#ifdef DEBUGGING
        BOOST_LOG_TRIVIAL(debug) << "Thread pool synced";
#endif // DEBUGGING
    }

    RNG * m_rng;
    sequence_bias_parameter< double > m_seq_bias;
    recombination_rate_parameter< double > m_recomb_rate;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_BATCH_CROSSOVER_MT_GENERATOR_HPP_
