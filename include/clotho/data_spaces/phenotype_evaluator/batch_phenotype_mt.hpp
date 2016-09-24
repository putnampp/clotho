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
#ifndef CLOTHO_BATCH_PHENOTYPE_MT_HPP_
#define CLOTHO_BATCH_PHENOTYPE_MT_HPP_

#include "clotho/data_spaces/phenotype_evaluator/phenotype_details.hpp"
#include "clotho/data_spaces/phenotype_evaluator/batch_phenotype_task.hpp"
#include "clotho/data_spaces/phenotype_evaluator/phenotype_accumulator.hpp"

namespace clotho {
namespace genetics {

template < class GeneticSpaceType >
class BatchPhenotypeMT : public phenotype_details< GeneticSpaceType > {
public:
    typedef phenotype_details< GeneticSpaceType >           base_type;
    typedef GeneticSpaceType                                genetic_space_type;

    typedef batch_phenotype_task< genetic_space_type >      task_type;

    BatchPhenotypeMT( boost::property_tree::ptree & config ) :
        base_type( config ) 
    {}

    template < class PoolType >
    void operator()( genetic_space_type * pop, PoolType & pool ) {
        this->resize( pop->getSequenceSpace().row_count(), this->m_trait_count );

        if( ! pop->getAlleleSpace().isAllNeutral() ) {
            generate( pop, pool );
            this->reduce_phenotypes(pop);
        }
    }

    virtual ~BatchPhenotypeMT() { }

protected:

    template < class PoolType >
    void generate( genetic_space_type * pop, PoolType & pool ) {
        unsigned int tc = pool.pool_size() + 1; // add one for master thread

        unsigned int batch_size = this->m_seq_count / tc;
        batch_size += ((this->m_seq_count % tc > 0) ? 1 : 0);

        unsigned int buffer_offset = batch_size * this->m_trait_count;

        typename base_type::phenotype_type * tmp = this->getPhenotypes();

        unsigned int i = 0;
        while( i + batch_size < this->m_seq_count ) {
            BOOST_LOG_TRIVIAL(info) << "Phenotype task (worker): [" << i << ", " << i + batch_size << "); Trait_count: " << this->m_trait_count;

            task_type t(pop, i, i + batch_size, this->m_trait_count, tmp );
            pool.post( t );

            i += batch_size;
            tmp += buffer_offset;
        }

        if( i < this->m_seq_count ) {
            BOOST_LOG_TRIVIAL(info) << "Phenotype task (master): [" << i << ", " << this->m_seq_count << "); Trait_count: " << this->m_trait_count;
            task_type t(pop, i, this->m_seq_count, this->m_trait_count, tmp );
            t();
        }

        pool.sync();
    }

};

}   // namespace genetics
}   // namespace clotho

namespace clotho {
namespace utility {

template < class GeneticSpaceType >
struct state_getter< clotho::genetics::BatchPhenotypeMT< GeneticSpaceType > >  {
    typedef clotho::genetics::BatchPhenotypeMT< GeneticSpaceType > object_type;

    void operator()( boost::property_tree::ptree & s, object_type & obj ) {
        typedef typename object_type::phenotype_type phenotype_type;

        size_t i = 0;
        phenotype_type  * tmp = obj.getPhenotypes();
        while( i < obj.sequence_count() ) {
            boost::property_tree::ptree p;
            clotho::utility::add_value_array( p, tmp, tmp + obj.trait_count() );

            s.push_back( std::make_pair( "", p ) );

            ++i;

            tmp += obj.trait_count();
        }
    }
};


}   // namespace utility
}   // namespace clotho

#endif  // CLOTHO_BATCH_PHENOTYPE_MT_HPP_

