#ifndef ALLELE_DISTRIBUTION_VARIABLE_HPP_
#define ALLELE_DISTRIBUTION_VARIABLE_HPP_

#include "clotho/genetics/allele_distribution_def.hpp"

#include "clotho/powerset/variable_subset.hpp"

template < class E, class B, class BM, class EK >
struct allele_distribution< clotho::powersets::variable_subset< E, B, BM, EK > > {

    typedef clotho::powersets::variable_subset< E, B, BM, EK > sequence_type;

    void update( sequence_type & seq );

};

#endif  // ALLELE_DISTRIBUTION_VARIABLE_HPP_
