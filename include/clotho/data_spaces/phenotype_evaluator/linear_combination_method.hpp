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
#ifndef CLOTHO_LINEAR_COMBINATION_METHOD_HPP_
#define CLOTHO_LINEAR_COMBINATION_METHOD_HPP_

#include "clotho/data_spaces/accumulator_helper_of.hpp"

namespace clotho {
namespace genetics {

template < class TraitAccumType, class ResultType >
struct linear_combination {

    typedef TraitAccumType                                  trait_accum_type;
    typedef typename trait_accum_type::trait_vector_type    trait_type;
    typedef ResultType                                      result_type;

    result_type operator()( trait_type & s0, trait_type & s1 ) {
        typedef typename trait_type::iterator iterator;

        iterator b = s0.begin(), e = s0.end();
        result_type res = 0.0;

        while( b != e ) {
            res += *b++;
        }

        b = s1.begin(); e = s1.end();

        while( b != e ) {
            res += *b++;
        }

        return res;
    }
};

template < class TraitAccumType, class ResultType >
struct accumulator_helper_of< linear_combination< TraitAccumType, ResultType > > {
    typedef TraitAccumType          type;
    typedef ResultType              result_type;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_LINEAR_COMBINATION_METHOD_HPP_
