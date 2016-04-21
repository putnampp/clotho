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
#ifndef CLOTHO_QTL_ALLELE_HPP_
#define CLOTHO_QTL_ALLELE_HPP_

#include "clotho/data_spaces/allele_space/neutral_allele.hpp"
#include "clotho/data_spaces/growable2D.hpp"
#include "clotho/data_spaces/trait_type_of.hpp"

namespace clotho {
namespace genetics {

template < class PositionType, class WeightType >
class qtl_allele_vectorized : public neutral_allele_vectorized< PositionType >, public growable2D {
public:
    typedef neutral_allele_vectorized< PositionType > base_type;

    typedef WeightType                  weight_type;
    typedef std::vector< WeightType >   trait_type;

    class weight_iterator {
    public:
        weight_iterator( weight_type * ptr, size_t size ) :
            m_ptr( ptr )
            , m_size( size ) {}

        weight_iterator( const weight_iterator & other ) :
            m_ptr( other.m_ptr )
            , m_size( other.m_size )
        {}

        bool hasNext() const {
            return (m_size > 0);
        }
        
        weight_type next() {
            assert( m_size > 0 );
            weight_type w = *m_ptr++;
            --m_size;

            return w;
        }

        virtual ~weight_iterator() {}
    protected:
        weight_type * m_ptr;
        size_t      m_size;
    };

    typedef weight_iterator     trait_iterator;

    qtl_allele_vectorized( size_t rows = 1, size_t columns = 1 ) :
        m_weights(NULL)
        , m_rows(0)
        , m_columns(0)
        , m_size(0)
    {
        this->resize(rows, columns);
    }
 
    trait_iterator getTraitIterator( size_t allele_idx ) const {
        assert( 0 <= allele_idx && allele_idx < m_rows );
        return trait_iterator( m_weights + allele_idx * m_columns, m_columns );
    }

    void updateTraitWeight( size_t allele_idx, size_t trait_idx, weight_type w ) {
        assert( 0 <= allele_idx && allele_idx < m_rows );
        assert( 0 <= trait_idx && trait_idx < m_columns );

        m_weights[ allele_idx * m_columns + trait_idx ] = w;
    }

    virtual size_t grow( size_t alleles ) {
        this->resize( alleles );
        return m_rows;
    }

    virtual size_t grow( size_t alleles, size_t traits ) {
        this->resize( alleles, traits );

        return m_rows;
    }

    trait_type  makeEmptyTrait() {
        return trait_type( m_columns, 0.0);
    }


    size_t allele_count() const {
        return m_rows;
    }

    size_t trait_count() const {
        return m_columns;
    }

    size_t allocated_size() const {
        return m_size;
    }

    virtual ~qtl_allele_vectorized() {
        if( m_weights != NULL ) {
            delete [] m_weights;
        }
    }

protected:

    virtual void resize( size_t rows, size_t columns ) {
        base_type::resize( rows );

        rows = base_type::size();

        size_t new_size = rows * columns;

        assert( 0 < new_size );

        if( m_size < new_size ) {
            if( m_weights != NULL ) {
                delete [] m_weights;
            }

            m_weights = new weight_type[ new_size ];

            m_size = new_size;
        }

        m_columns = columns;
        m_rows = rows;
    }

    virtual void resize( size_t rows ) {
        this->resize( rows, m_columns );
    }

    weight_type * m_weights;
    size_t      m_rows, m_columns, m_size;
};

template < class PositionType, class WeightType >
struct trait_type_of< qtl_allele_vectorized< PositionType, WeightType > > {
    typedef qtl_allele_vectorized< PositionType, WeightType > allele_type;
    typedef typename allele_type::trait_type type;

    static type makeEmptyTrait( size_t N ) {
        type r = type( N, 0.0);
        return r;
    }

    static void resetTrait( type & trait ) {
        size_t i = 0;
        while( i < trait.size() ) {
            trait[i++] = 0.0;
        }
    }

};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_QTL_ALLELE_HPP_