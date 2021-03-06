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
#ifndef INDIVIDUAL_TRAIT_ACCUMULATOR_HPP_
#define INDIVIDUAL_TRAIT_ACCUMULATOR_HPP_

#include "sequence_helper.hpp"
#include "clotho/genetics/sequence_trait_accumulator.hpp"

template < class IndividualType >
class individual_trait_accumulator {
public:
    typedef IndividualType                                  individual_type;
    typedef sequence_helper< individual_type >              sequence_helper_type;
    typedef typename sequence_helper_type::sequence_type    sequence_type;

    typedef sequence_trait_accumulator< sequence_type >     sequence_accumulator_type;
    typedef typename sequence_accumulator_type::result_type sequence_result_type;
    typedef typename sequence_accumulator_type::result_type result_type;

    typedef typename sequence_helper_type::iterator_helper_type   iterator_helper_type;
    typedef typename iterator_helper_type::type             iterator;

    individual_trait_accumulator( result_type & res ) : m_acc( res ) {}

    inline void operator()( sequence_result_type & sres ) {
        m_acc( sres );
    }

    result_type operator()( individual_type & ind ) {
        iterator first = iterator_helper_type::make_first(ind), last = iterator_helper_type::make_last(ind);

        result_type res;

        sequence_accumulator_type acc( res );
        std::for_each( first, last, acc );

        return res;
    }

protected:
    sequence_accumulator_type   m_acc;
};

#endif  // INDIVIDUAL_TRAIT_ACCUMULATOR_HPP_
