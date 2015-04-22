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
#ifndef FITNESS_CORE_HPP_
#define FITNESS_CORE_HPP_

#include "clotho/fitness/policies/fitness_policies.hpp"
#include "clotho/fitness/support/allele_set_of.hpp"

namespace clotho {
namespace fitness {
namespace detail {


template < class Result, class Individual, class SelectionTag >
class fitness_core {
public:
    typedef Result result_type;

    typedef clotho::fitness::support::allele_set_of< Individual >::type           set_type;
    typedef clotho::fitness::support::allele_set_of< Individual >::subset_type    subset_type;
    typedef clotho::fitness::support::allele_set_of< Individual >::element_type   allele_type;

    typedef clotho::fitness::policies::fitness_policy< result_type, allele_type, SelectionTag > fitness_policy_type;

    typedef clotho::fitness::detail::eval_method< result_type, subset_type, fitness_policy_type > evaluation_type;

};

}   // namespace detail
}   // namespace fitness
}   // namespace clotho

#endif  // FITNESS_CORE_HPP_
