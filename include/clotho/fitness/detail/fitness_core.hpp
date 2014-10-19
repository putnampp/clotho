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
