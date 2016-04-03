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
#ifndef EFFECT_SIZE_MATRIX_HPP_
#define EFFECT_SIZE_MATRIX_HPP_

template < class EffSizeType >
class EffectSizeMatrix {
public:
    typedef EffSizeType value_type;

    value_type *    m_eff_size_mat;

    unsigned int    m_allele_count, m_trait_count;
    unsigned int    m_size;

    EffectSizeMatrix( unsigned int allele_count, unsigned int trait_count ) :
        m_eff_size_mat( NULL )
        , m_allele_count( 0 )
        , m_trait_count( 0 )
        , m_size( 0 )
    {
        resize( allele_count, trait_count);
    }

    template< class AlleleType >
    bool update( const AlleleType & alleles ) {
        bool all_neutral = true;
        unsigned int i = 0;

        unsigned int a_count = alleles.size();

        typename AlleleType::cvariable_iterator it = alleles.variable_begin();

        unsigned int t_count = it->trait_count();

        resize( a_count, t_count );

        while( it != alleles.variable_end() ) {
            all_neutral = all_neutral && it->isNeutral();
            for( typename AlleleType::value_type::weight_citerator wit = it->begin(); wit != it->end(); wit++) {
                m_eff_size_mat[ i++ ] = *wit;
            }
            it++;
        }

        return all_neutral;
    }

    virtual ~EffectSizeMatrix() {
        if( m_eff_size_mat != NULL ) {
            delete m_eff_size_mat;
        }
    }

protected:
    void resize( unsigned int allele_count, unsigned int trait_count ) {
        unsigned int new_size = allele_count * trait_count;

        if( new_size > m_size ) {
            // growing
            //
            if( m_eff_size_mat != NULL ) {
                delete m_eff_size_mat;
            }
            
            m_eff_size_mat = new value_type[ new_size ];
            m_size = new_size;
        }
        m_allele_count = allele_count;
        m_trait_count = trait_count;
    }
};
#endif  // EFFECT_SIZE_MATRIX_HPP_
