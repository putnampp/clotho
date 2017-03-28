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
#ifndef CLOTHO_TRAIT_SPACE_VECTOR_HPP_
#define CLOTHO_TRAIT_SPACE_VECTOR_HPP_

#include "clotho/data_spaces/phenotype_evaluator/trait_count_parameter.hpp"
#include <vector>

#include "clotho/utility/state_object.hpp"
#include "clotho/utility/log_helper.hpp"

namespace clotho {
namespace genetics {

template < class WeightType >
class trait_space_vector {
public:
    typedef WeightType       weight_type;

    typedef std::vector< weight_type > vector_type;
    typedef typename vector_type::iterator iterator;
    typedef typename vector_type::const_iterator const_iterator;

    trait_space_vector( boost::property_tree::ptree & config ) :
        m_trait_count(0)
    {
        trait_count_parameter tc_param( config );
        m_trait_count = tc_param.m_trait_count;
    }

    unsigned int trait_count() const {
        return m_trait_count;
    }

    unsigned int weight_count() const {
        return m_weights.size();
    }

    unsigned int allele_count() const {
        return m_weights.size() / m_trait_count;
    }

    weight_type operator[]( unsigned int offset ) const {
#ifdef DEBUGGING
        assert( offset < m_weights.size() );
#endif  // DEBUGGING
        return m_weights[ offset ];
    }

    weight_type getWeight( unsigned int row, unsigned int col ) const {
        row *= m_trait_count;
        row += col;

#ifdef DEBUGGING
        assert( row < m_weights.size() );
#endif  // DEBUGGING
        return m_weights[ row ];
    }

//    void setWeight( weight_type w, unsigned int row, unsigned int col ) {
//#ifdef DEBUGGING
//        assert( col < m_trait_count );
//#endif  // DEBUGGING
//
//        unsigned int _offset = row * m_trait_count;
//
//        if( _offset < m_weights.size() ) {
//            m_weights[ _offset + col ] = w;
//        } else {
//            do {
//                for( unsigned int i = 0; i < m_trait_count; ++i ) {
//                    m_weights.push_back( w );
//                }
//            } while( m_weights.size() <= _offset );
//        }
//
//    }
//
    void setWeight( weight_type w, unsigned int row, unsigned int col ) {
        assert( col < m_trait_count );
        unsigned int _offset = row * m_trait_count;

//        while( m_weights.size() <= _offset ) grow();

        m_weights[ _offset + col ] = w;
    }

    void append( iterator first, iterator last ) {
        while( first != last ) {
            m_weights.push_back( *first++ );
        }
    }

    void grow() {
        for( unsigned int i = 0; i < m_trait_count; ++i ) {
            m_weights.push_back( 0.0 );
        }
    }

    iterator begin( unsigned int r ) {
        r *= m_trait_count;

#ifdef DEBUGGING
        assert( r < m_weights.size() );
#endif  // DEBUGGING
        return m_weights.begin() + r;
    }

    iterator end( unsigned int r ) {
        ++r;
        r *= m_trait_count;

        if( r < m_weights.size() ) {
            return m_weights.begin() + r;
        }
#ifdef DEBUGGING
         else if ( r == m_weights.size() ) {
            return m_weights.end();
        } else {
            assert( false );
        }
#else
        else {
            return m_weights.end();
        }
#endif  // DEBUGGING
    }

    const_iterator begin( unsigned int r ) const {
        r *= m_trait_count;
#ifdef DEBUGGING
        assert( r < m_weights.size() );
#endif  // DEBUGGING

        return m_weights.begin() + r;
    }

    const_iterator end( unsigned int r ) const {
        ++r;
        r *= m_trait_count;

        if( r < m_weights.size() ) {
            return m_weights.begin() + r;
        }
#ifdef DEBUGGING
         else if ( r == m_weights.size() ) {
            return m_weights.end();
        } else {
            // if this was encountered, then there is some bad math upstream of the call.
            // fix it
            assert( false );
        }
#else
        else {
            return m_weights.end();
        }
#endif  // DEBUGGING
    }
    
    virtual ~trait_space_vector() {}

protected:
    unsigned int m_trait_count;

    vector_type m_weights;
};

}   // namespace genetics
}   // namespace clotho

namespace clotho {
namespace utility {

template < class WeightType >
struct state_getter< clotho::genetics::trait_space_vector< WeightType > > {
    typedef clotho::genetics::trait_space_vector< WeightType > object_type;

    void operator()( boost::property_tree::ptree & s, object_type & obj ) {
        typedef typename object_type::iterator    iterator;

        size_t all_count = obj.allele_count();
        size_t i = 0;
        while( i < all_count ) {

            boost::property_tree::ptree t;
            iterator first = obj.begin( i ), last = obj.end( i );

            clotho::utility::add_value_array( t, first, last );

            s.push_back( std::make_pair("",  t ) );
            ++i;
        }

    }
};

}   // namespace utility
}   // namespace clotho

#endif  // CLOTHO_TRAIT_SPACE_VECTOR_HPP_

