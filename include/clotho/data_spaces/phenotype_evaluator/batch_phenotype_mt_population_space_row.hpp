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
#ifndef CLOTHO_BATCH_PHENOTYPE_MT_POPULATION_SPACE_ROW_HPP_
#define CLOTHO_BATCH_PHENOTYPE_MT_POPULATION_SPACE_ROW_HPP_

#include "clotho/data_spaces/population_space/population_space_row.hpp"

#include "clotho/data_spaces/phenotype_evaluator/batch_phenotype_task_population_space_row.hpp"
#include "clotho/data_spaces/phenotype_evaluator/individual_reducer_tasks.hpp"

namespace clotho {
namespace genetics {

template < class BlockType, class WeightType, class TraitSpaceType >
class BatchPhenotypeMT< population_space_row< BlockType, WeightType >, TraitSpaceType > {
public:

    typedef BatchPhenotypeMT< population_space_row< BlockType, WeightType >, TraitSpaceType > self_type;
    typedef population_space_row< BlockType, WeightType >                      space_type;
    typedef TraitSpaceType                                                          trait_space_type;

    typedef typename space_type::weight_type                                        phenotype_type;
    typedef typename space_type::weight_vector                                      phenotype_vector;

    typedef typename space_type::weight_vector                                      phenotype_space_type;

    typedef typename phenotype_space_type::iterator                                 phenotype_iterator;
    typedef typename phenotype_space_type::const_iterator                           const_phenotype_iterator;

    typedef batch_phenotype_task< space_type, trait_space_type >                    updater_task;
    typedef individual_reducer_task< space_type, self_type >                        reducer_type;

    BatchPhenotypeMT() : 
        m_trait_count(0) 
    {}

    void constant_phenotype( space_type * pop, trait_space_type * traits ) {
        m_trait_count = pop->getTraitCount();
        unsigned int i = 0, pheno_size = pop->getTraitCount() * pop->getIndividualCount();

        while( i < m_phenos.size() ) {
            m_phenos[ i++ ] = 1.0;
        }

        while( m_phenos.size() < pheno_size ) {
            m_phenos.push_back(1.0);
        }
    }

    template < class PoolType >
    void operator()( space_type * pop, trait_space_type * traits, PoolType & pool ) {
        m_trait_count = pop->getTraitCount();

        const unsigned int N = pop->getTraitCount() * pop->getIndividualCount();
        while( m_phenos.size() < N ) {
            m_phenos.push_back(0.0);
        }

        updateHaploidGenomeWeights( pop, traits, pool );
        updateIndividualPhenotypes( pop, traits, pool );
    }

    phenotype_iterator begin() {
        return m_phenos.begin();
    }

    phenotype_iterator end() {
        return m_phenos.end();
    }

    const_phenotype_iterator begin() const {
        return m_phenos.begin();
    }

    const_phenotype_iterator end() const {
        return m_phenos.end();
    }

    phenotype_iterator begin_individual_phenotype( unsigned int i ) {
        return m_phenos.begin() + i * m_trait_count;
    }

    phenotype_iterator end_individual_phenotype( unsigned int i ) {
        return m_phenos.begin() + (i + 1) * m_trait_count;
    }

    unsigned int individual_count() const {
        return m_phenos.size();
    }

    unsigned int trait_count() const {
        return m_trait_count;
    }

    virtual ~BatchPhenotypeMT() {}

protected:

    template < class PoolType >
    void updateHaploidGenomeWeights( space_type * pop, trait_space_type * traits, PoolType & pool ) {
        unsigned int tc = pool.pool_size() + 1;

        const unsigned int M = pop->haploid_genome_count();
        unsigned int batch_size = M / tc + ((M % tc > 0) ? 1 : 0);

        unsigned int i = 0;
        while( i + batch_size < M ) {

            updater_task t( pop, traits, i, i + batch_size );
            t();

            i += batch_size;
        }

        if( i < M ) {
            updater_task t( pop, traits, i, M );
            t();
        }

        pool.sync();
    }

    template < class PoolType >
    void updateIndividualPhenotypes( space_type * pop, trait_space_type * traits, PoolType & pool ) {
        unsigned int tc = pool.pool_size() + 1;

        const unsigned int M = pop->getIndividualCount();
        unsigned int batch_size = M / tc + ((M % tc > 0) ? 1 : 0);

        unsigned int i = 0;
        while( i + batch_size < M ) {

            reducer_type t( pop, i, i + batch_size, this );
            t();

            i += batch_size;
        }

        if( i < M ) {
            reducer_type t( pop, i, M, this );
            t();
        }

        pool.sync();
    }

    unsigned int            m_trait_count;
    phenotype_space_type    m_phenos;
};

}   // namespace genetics
}   // namespace clotho

namespace clotho {
namespace utility {

template < class BlockType, class WeightType, class TraitSpaceType >
struct state_getter< clotho::genetics::BatchPhenotypeMT< clotho::genetics::population_space_row<BlockType, WeightType>, TraitSpaceType > >  {
    typedef clotho::genetics::BatchPhenotypeMT< clotho::genetics::population_space_row<BlockType, WeightType>, TraitSpaceType > object_type;

    void operator()( boost::property_tree::ptree & s, object_type & obj ) {
//        typedef typename object_type::phenotype_type phenotype_type;
//
//        size_t i = 0;
//        phenotype_type  * tmp = obj.getPhenotypes();
//        while( i < obj.sequence_count() ) {
//            boost::property_tree::ptree p;
//            clotho::utility::add_value_array( p, tmp, tmp + obj.trait_count() );
//
//            s.push_back( std::make_pair( "", p ) );
//
//            ++i;
//
//            tmp += obj.trait_count();
//        }
    }
};


}   // namespace utility
}   // namespace clotho
#endif  // CLOTHO_BATCH_PHENOTYPE_MT_POPULATION_SPACE_ROW_HPP_
