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
#ifndef CLOTHO_PHENOTYPE_MT_HPP_
#define CLOTHO_PHENOTYPE_MT_HPP_

#include "clotho/data_spaces/phenotype_evaluator/phenotype_details.hpp"
#include "clotho/data_spaces/phenotype_evaluator/phenotype_task.hpp"
#include "clotho/data_spaces/phenotype_evaluator/phenotype_accumulator.hpp"

namespace clotho {
namespace genetics {

template < class SequenceSpaceType, class TraitSpaceType >
class PhenotypeMT : public phenotype_details< SequenceSpaceType, TraitSpaceType > {
public:

    typedef phenotype_details< SequenceSpaceType, TraitSpaceType >   base_type;
    typedef SequenceSpaceType                                       sequence_space_type;
    typedef TraitSpaceType                                          trait_space_type;

    typedef typename base_type::weight_type                         weight_type;
    typedef typename sequence_space_type::block_type                block_type;
    typedef typename sequence_space_type::sequence_vector           sequence_vector;

    typedef phenotype_accumulator< trait_space_type >               accumulator_type;
    typedef phenotype_task< block_type, accumulator_type >          task_type;

    PhenotypeMT( boost::property_tree:: ptree & config ) :
        base_type( config )
    {    }

    template < class PoolType >
    void operator()( sequence_space_type * pop, trait_space_type * traits, PoolType & pool ) {
        this->resize( pop->row_count(), traits->trait_count() );

        if( pool.pool_size() < 1 ) {
            generate( pop, traits );
        } else {
            generate( pop, traits, pool );
        }
        this->reduce_phenotypes();
    }

    virtual ~PhenotypeMT() { }

protected:

    void generate( sequence_space_type * pop, trait_space_type * traits ) {
        typename base_type::phenotype_type * tmp = this->getPhenotypes();
        for( unsigned int i = 0; i < this->m_seq_count; ++i ) {
            accumulator_type acc( traits, tmp );

            sequence_vector seq = pop->getSequence( i );

            task_type t( seq.first, seq.second, acc );

            t();

            tmp += traits->trait_count();
        }
    }

    template < class PoolType >
    void generate( sequence_space_type * pop, trait_space_type * traits, PoolType & pool ) {
        typename base_type::phenotype_type * tmp = this->getPhenotypes();
        for( unsigned int i = 0; i < this->m_seq_count; ++i ) {
            accumulator_type acc( traits, tmp );

            sequence_vector seq = pop->getSequence( i );

            task_type t(seq.first, seq.second, acc );
            pool.post( t );

            tmp += traits->trait_count();
        }

        pool.sync();
    }
};

}   // namespace genetics
}   // namespace clotho

namespace clotho {
namespace utility {

template < class SequenceSpaceType, class TraitSpaceType >
struct state_getter< clotho::genetics::PhenotypeMT< SequenceSpaceType, TraitSpaceType > >  {
    typedef clotho::genetics::PhenotypeMT< SequenceSpaceType, TraitSpaceType > object_type;

    void operator()( boost::property_tree::ptree & s, object_type & obj ) {
        typedef typename object_type::weight_type phenotype_type;

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
#endif  // CLOTHO_PHENOTYPE_MT_HPP_
