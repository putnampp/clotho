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

template < class SequenceSpaceType, class TraitSpaceType >
class BatchPhenotypeMT : public phenotype_details< SequenceSpaceType, TraitSpaceType > {
public:
    typedef phenotype_details< SequenceSpaceType, TraitSpaceType >              base_type;
    typedef SequenceSpaceType                                                   sequence_space_type;
    typedef TraitSpaceType                                                      trait_space_type;

    typedef batch_phenotype_task< sequence_space_type, trait_space_type >       task_type;

    BatchPhenotypeMT( ) {}

    template < class PoolType >
    void operator()( sequence_space_type * pop, trait_space_type * trait_space, PoolType & pool ) {
        this->resize( pop->row_count(), trait_space->trait_count() );

        evaluate( pop, trait_space, pool );
        this->reduce_phenotypes();
    }

    virtual ~BatchPhenotypeMT() { }

protected:

    template < class PoolType >
    void evaluate( sequence_space_type * pop, trait_space_type * traits, PoolType & pool ) {
        unsigned int tc = pool.pool_size() + 1; // add one for master thread

        unsigned int batch_size = this->m_seq_count / tc;
        batch_size += ((this->m_seq_count % tc > 0) ? 1 : 0);

        unsigned int buffer_offset = batch_size * traits->trait_count();

        typename base_type::phenotype_type * tmp = this->getPhenotypes();

        unsigned int i = 0;
        while( i + batch_size < this->m_seq_count ) {
            BOOST_LOG_TRIVIAL(info) << "Phenotype task (worker): [" << i << ", " << i + batch_size << "); Trait_count: " << this->m_trait_count;

            task_type t(pop, traits, i, i + batch_size, tmp );
            pool.post( t );

            i += batch_size;
            tmp += buffer_offset;
        }

        if( i < this->m_seq_count ) {
            BOOST_LOG_TRIVIAL(info) << "Phenotype task (master): [" << i << ", " << this->m_seq_count << "); Trait_count: " << this->m_trait_count;
            task_type t(pop, traits, i, this->m_seq_count, tmp );
            t();
        }

        pool.sync();
    }

};

}   // namespace genetics
}   // namespace clotho

namespace clotho {
namespace utility {

template < class SequenceSpaceType, class TraitSpaceType >
struct state_getter< clotho::genetics::BatchPhenotypeMT< SequenceSpaceType, TraitSpaceType > >  {
    typedef clotho::genetics::BatchPhenotypeMT< SequenceSpaceType, TraitSpaceType > object_type;

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

