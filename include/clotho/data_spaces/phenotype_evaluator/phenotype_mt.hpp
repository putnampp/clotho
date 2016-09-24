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

template < class GeneticSpaceType >
class PhenotypeMT : public phenotype_details< GeneticSpaceType > {
public:

    typedef phenotype_details< GeneticSpaceType >   base_type;
    typedef GeneticSpaceType                        genetic_space_type;

    typedef typename base_type::weight_type                         weight_type;
    typedef typename genetic_space_type::block_type                 block_type;
    typedef typename genetic_space_type::association_type::sequence_vector  sequence_vector;

    typedef phenotype_accumulator< weight_type >                    accumulator_type;
    typedef phenotype_task< block_type, accumulator_type >          task_type;

    PhenotypeMT( boost::property_tree:: ptree & config ) :
        base_type( config )
    {    }

    template < class PoolType >
    void operator()( genetic_space_type * pop, PoolType & pool ) {
        this->resize( pop->getSequenceSpace().row_count(), this->m_trait_count );

        if( !pop->getAlleleSpace().isAllNeutral() ) {
            if( pool.pool_size() < 1 ) {
                generate( pop );
            } else {
                generate( pop, pool );
            }
            this->reduce_phenotypes(pop);
        }
    }

    virtual ~PhenotypeMT() { }

protected:

    void generate( genetic_space_type * pop ) {
        typename base_type::phenotype_type * tmp = this->getPhenotypes();
        for( unsigned int i = 0; i < this->m_seq_count; ++i ) {
            accumulator_type acc( pop->getAlleleSpace().getWeights(), pop->getAlleleSpace().allele_count(), this->m_trait_count, tmp );

            sequence_vector seq = pop->getSequenceSpace().getSequence( i );

            task_type t( seq.first, seq.second, acc );

            t();

            tmp += this->m_trait_count;
        }
    }

    template < class PoolType >
    void generate( genetic_space_type * pop, PoolType & pool ) {
        typename base_type::phenotype_type * tmp = this->getPhenotypes();
        for( unsigned int i = 0; i < this->m_seq_count; ++i ) {
            accumulator_type acc( pop->getAlleleSpace().getWeights(), pop->getAlleleSpace().allele_count(), this->m_trait_count, tmp );

            sequence_vector seq = pop->getSequenceSpace().getSequence( i );

            task_type t(seq.first, seq.second, acc );
            pool.post( t );

            tmp += this->m_trait_count;
        }

        pool.sync();
    }
};

}   // namespace genetics
}   // namespace clotho

namespace clotho {
namespace utility {

template < class GeneticSpaceType >
struct state_getter< clotho::genetics::PhenotypeMT< GeneticSpaceType > >  {
    typedef clotho::genetics::PhenotypeMT< GeneticSpaceType > object_type;

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
