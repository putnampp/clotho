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
#ifndef CLOTHO_PHENOTYPE_DETAILS_HPP_
#define CLOTHO_PHENOTYPE_DETAILS_HPP_

#include "clotho/data_spaces/phenotype_evaluator/pairwise_phenotype_reducer.hpp"

namespace clotho {
namespace genetics {

template < class SequenceSpaceType, class TraitSpaceType >
class phenotype_details {
public:
    typedef SequenceSpaceType                                       sequence_space_type;
    typedef TraitSpaceType                                          trait_space_type;

    typedef typename trait_space_type::weight_type                  weight_type;
    typedef weight_type                                             phenotype_type;

    phenotype_details( ) :
        m_phenos( NULL )
        , m_seq_count( 0 )
        , m_trait_count( 0 )
        , m_alloc_size( 0 ) 
    { }

    void constant_phenotype( sequence_space_type * pop, trait_space_type * trait_space ) {
        resize( pop->row_count(), trait_space->trait_count() );

        for( unsigned int i = 0; i < m_alloc_size; ++i ) {
            m_phenos[ i ] = 1.0;
        }
    }

    phenotype_type *   getPhenotypes() const {
        return m_phenos;
    }

    unsigned int individual_count() const {
        return m_seq_count / 2;
    }

    unsigned int sequence_count() const {
        return m_seq_count;
    }
    
    unsigned int trait_count() const {
        return m_trait_count;
    }

    virtual ~phenotype_details() {
        if( m_phenos != NULL ) {
            delete [] m_phenos;
        }
    }

protected:

    void reduce_phenotypes(  ) {
        pairwise_phenotype_reducer red;

        unsigned int i = 0;
        phenotype_type * res = this->m_phenos;
        phenotype_type * tmp = this->m_phenos;

        while( i < this->m_seq_count ) {
            red( res, tmp, tmp + this->m_trait_count, this->m_trait_count );

            res += this->m_trait_count;
            tmp += 2 * this->m_trait_count;
            i += 2;  
        }
    }

    void resize( unsigned int seq_count, unsigned int trait_count ) {

        unsigned int s = seq_count * trait_count;
        if( s > m_alloc_size ) {
            if( m_phenos != NULL ) {
                delete [] m_phenos;
            }

            m_phenos = new weight_type[ s ];

            m_alloc_size = s;
        }

        memset( m_phenos, 0, m_alloc_size * sizeof(weight_type) );

        m_seq_count = seq_count;
        m_trait_count = trait_count;
    }

    weight_type     * m_phenos;
    unsigned int    m_seq_count, m_trait_count, m_alloc_size;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_PHENOTYPE_DETAILS_HPP_

