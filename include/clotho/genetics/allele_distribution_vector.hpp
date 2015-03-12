#ifndef ALLELE_DISTRIBUTION_VECTOR_HPP_
#define ALLELE_DISTRIBUTION_VECTOR_HPP_

#include "clotho/genetics/allele_distribution_def.hpp"

#include "clotho/powerset/vector_subset.hpp"

template < class E, class B, class BM, class EK >
struct allele_distribution < clotho::powersets::vector_subset< E, B, BM, EK > > {
    typedef clotho::powersets::vector_subset< E, B, BM, EK > sequence_type;

    void update( sequence_type & seq );

};
#endif  //ALLELE_DISTRIBUTION_VECTOR_HPP_
