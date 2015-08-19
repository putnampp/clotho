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
#ifndef HASH_CLASSIFIER_HPP_
#define HASH_CLASSIFIER_HPP_

#include <vector>
#include <cstring>
#include <functional>
#include "clotho/classifiers/helpers/normalize.hpp"

namespace clotho {
namespace classifiers {

template < class Element, unsigned int BinCount = 256 >
class hash_classifier;

template < class Element, unsigned int BinCount >
class hash_classifier {
public:
    typedef hash_classifier< Element, BinCount > self_type;
    typedef Element element_type;
    typedef size_t  result_type;

    typedef std::hash< self_type > hash_type;

    static const unsigned int BINS = BinCount;

    struct param {
        typedef unsigned int bin_type;

        bin_type bin_psum[ BINS + 1 ];  // prefix sum of class_bounds per bin
        std::vector< element_type > class_bounds;

        param( ) {
            memset( bin_psum, 0, sizeof( bin_type ) * (BINS + 1) );
        }

        param( const param & p ) : class_bounds( p.class_bounds ) {
            memcpy( bin_psum, p.bin_psum, sizeof( bin_type ) * (BINS + 1) );
        }
    };

    typedef param param_type;

    hash_classifier() {}

    hash_classifier( const param_type & p ) : m_param( p ) {}

    hash_classifier( const self_type & rhs ) : m_param( rhs.m_param ) {}

    param_type * param() {
        return &m_param;
    }

    /// returns the class index
    inline result_type operator()( const element_type & elem ) const {
//        size_t bin_idx = to_key< element_type, double >::convert( elem ) * BINS;
        size_t bin_idx = m_lookup( elem );

        unsigned int c_start = m_param.bin_psum[ bin_idx++ ];
        unsigned int c_end = m_param.bin_psum[ bin_idx ];

        while( c_start < c_end ) {
            if( !(m_param.class_bounds[ c_start ] < elem) ) {
                break;
            }
            ++c_start;
        }

        return c_start;
    }

    template < class Iter >
    result_type operator()( Iter first, size_t offset ) {
        return operator()( *(first + offset));
    }

    virtual ~hash_classifier() {}
protected:
    param_type  m_param;
    hash_type m_lookup;
};

}   // namespace classifiers
}   // namespace clotho


namespace std {

template < class Element, unsigned int Bins >
struct hash< clotho::classifiers::hash_classifier< Element, Bins > > {
    typedef Element     argument_type;
    typedef std::size_t result_type;

    result_type operator()( argument_type const& elem ) {
        result_type res = ( normalize_arg( elem ) * Bins);
        return res;
    }
};

}   // namespace std

#endif  // HASH_CLASSIFIER_HPP_
