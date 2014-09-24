#ifndef FITNESS_GUARD_HPP_
#define FITNESS_GUARD_HPP_

//#include "clotho/fitness/fitness_bitset.hpp"
//

#include "clotho/fitness/detail/fitness_core.hpp"

template < class Result, class Individual, class SelectionTag >
class fitness :
public clotho::fitness::detail::fitness_core< Result, Individual, SelectionTag >,
public clotho::fitness::detail::fitness_core< Result, Individual, SelectionTag >::evaluation_type
{
public:
    
    typedef clotho::fitness::detail::fitness_core< Result, Individual, SelectionTag > core;
    typedef core::evaluation_type eval_type;

    typedef eval_type::hom_policy_type  hom_policy_type;
    typedef eval_type::het_policy_type  het_policy_type;

    fitness( const hom_policy_type & hom, const het_policy_type & het ) : eval_type( hom, het ) {}

};


#endif  // FITNESS_GUARD_HPP_
