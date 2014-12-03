#ifndef TRAIT_WEIGHT_HPP_
#define TRAIT_WEIGHT_HPP_

#include <vector>

extern const string TRAIT_BLOCK_K;

template < class ValueType = double >
struct trait_weight {
    typedef ValueType                   value_type;
    typedef std::vector< value_type >   vector_type;
};

#include "trait_weight_generator.hpp"

#endif  // TRAIT_WEIGHT_HPP_
