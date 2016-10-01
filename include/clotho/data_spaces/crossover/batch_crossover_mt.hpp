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

namespace clotho {
namespace genetics {

template < class RNG, class SequenceSpaceType, class AlleleSpaceType >
class BatchCrossoverMT {
public:
    typedef SequenceSpaceType       sequence_space_type;
    typedef AlleleSpaceType         allele_space_type;

//    typedef batch_crossover_task< RNG, SequenceSpaceType, AlleleSpaceType >     crossover_task_type;

    BatchCrossoverMT( RNG * rng, boost::property_tree::ptree & config ) :
        m_rng( rng )
        , m_seq_bias(config)
        , m_recomb_rate(config)
    { }

    template < class SelectionType, class PoolType >
    void operator()( const SelectionType & parents, sequence_space_type * parental, sequence_space_type * offspring, allele_space_type * alleles, PoolType & pool ) {
        offspring->clear();

        batch_generate( parents, parental, offspring, alleles, pool );
    }

    virtual ~BatchCrossoverMT() {}

protected:
    
    template < class SelectionType, class PoolType >
    void batch_generate( const SelectionType & parents, sequence_space_type * parental, sequence_space_type * offspring, allele_space_type * alleles, PoolType & pool ) {
#ifdef  DEBUGGING
        BOOST_LOG_TRIVIAL(debug) << "Launching crossover batch jobs";
#endif  // DEBUGGING
        typedef typename SelectionType::mate_pair_vector                            mate_pair_type;
        typedef batch_crossover_task< RNG, mate_pair_type, SequenceSpaceType, AlleleSpaceType >     crossover_task_type;
    
        size_t tc = pool.pool_size() + 1;   // + 1 for master thread

        size_t  batch_size = parents.size() / tc;
        batch_size += ((parents.size() % tc > 0)? 1 : 0);

        unsigned int off_idx = 0, j = 0;
        while( off_idx + batch_size < parents.size() ) {
            unsigned int off_end = off_idx + batch_size;

            BOOST_LOG_TRIVIAL(info) << "Batch Crossover: [" << off_idx << ", " << off_end << ");";

            crossover_task_type x( pool.getRNG(j++), parental, offspring, alleles, off_idx, parents.begin() + off_idx, parents.begin() + off_end, m_recomb_rate.m_rho, m_seq_bias.m_bias );
            pool.post( x );

            off_idx = off_end;
        }

        if( off_idx < parents.size() ){
            BOOST_LOG_TRIVIAL(info) << "Batch Crossover: [" << off_idx << ", " << parents.size() << ");";
            crossover_task_type t( m_rng, parental, offspring, alleles, off_idx, parents.begin() + off_idx, parents.end(), m_recomb_rate.m_rho, m_seq_bias.m_bias);
            t();
        }

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
