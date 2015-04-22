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
#ifndef INDIVIDUAL_GENOTYPER_HPP_
#define INDIVIDUAL_GENOTYPER_HPP_

#include <vector>
#include <algorithm>

#include "iterator_helper.hpp"

struct add {
    template < class ResultType >
    ResultType operator()( ResultType a, ResultType b ) {
        return a + b;
    }
};

template < class SetType, class ValueType = double, class OP = add >
class individual_genotyper {
public:
    typedef individual_genotyper< SetType, ResultType > self_type;
    typedef ValueType       value_type;
    typedef SetType         set_type;

    typedef typename set_type::iterator iterator;

    class result_type {
    public:
        typedef std::vector< value_type > values_type;
        typedef typename values_type::iterator iterator;
        typedef typename values_type::const_iterator citerator;

        friend class individual_genotyper< SetType, ValueType, OP >;

        result_type ( size_t n = 1 ) : m_values( n ) {}
        result_type ( const result_type & r ) : m_values( r.m_values ) {}

        iterator    begin() {
            return m_values.begin();
        }
        citerator   begin() const {
            return m_values.begin();
        }

        iterator    end() {
            return m_values.end();
        }
        citerator   end() const {
            return m_values.end();
        }

    protected:
        values_type m_values;
    };

    individual_genotyper( set_type & s, size_t n = 1 ) : m_set( &s ), m_nTraits(n) {}

    template < class Sequence >
    result_type operator()( Sequence & seq ) {

        typedef iterator_helper< Sequence >      allele_helper;
        typedef typename allele_helper::type     allele_iterator;

        allele_iterator first = allele_helper::make_first( seq ), last = allele_helper::make_last( seq );

        std::for_each( first, last, weight );
    }

protected:
    set_type    * m_set;
    size_t      m_nTraits;
};

#endif  // INDIVIDUAL_GENOTYPER_HPP_
