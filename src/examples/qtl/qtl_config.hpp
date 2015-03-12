#ifndef QTL_ESTIMATOR_CONFIG_HPP_
#define QTL_ESTIMATOR_CONFIG_HPP_

#include <boost/mpl/aux_/na.hpp>

#ifdef USE_VECTOR_SUBSET
#include "clotho/powerset/vector_subset.hpp"
#define SUBSETTYPE clotho::powersets::vector_subset

#include "clotho/recombination/inspect_methods.hpp"

#define RECOMBINE_INSPECT_METHOD clotho::recombine::inspection::tag::copy_matching_classify_mismatch
#define BIT_WALK_METHOD boost::mpl::na

#else // ! USE_VECTOR_SUBSET

#include "clotho/powerset/variable_subset.hpp"
#define SUBSETTYPE clotho::powersets::variable_subset

#if defined( USE_SCAN_AND_CLASSIFY )
#define BIT_WALK_METHOD clotho::recombine::walker::tag::scan_and_classify
#elif  defined( USE_SCAN_AND_CLASSIFY_SWITCH )
#define BIT_WALK_METHOD clotho::recombine::walker::tag::scan_and_classify_switch
#elif  defined( USE_INLINE_AND_CLASSIFY )
#define BIT_WALK_METHOD clotho::recombine::walker::tag::inline_classify
#elif   defined( USE_INLINE_DYNAMIC_AND_CLASSIFY )
#define BIT_WALK_METHOD clotho::recombine::walker::tag::inline_dynamic_classify
#else   // default is to perform classification inline
#define BIT_WALK_METHOD clotho::recombine::walker::tag::iterate_and_classify
#endif  // classification procedure

#if defined( USE_INSPECT_ALL )
#define RECOMBINE_INSPECT_METHOD clotho::recombine::inspection::tag::classify_all
#else
#define RECOMBINE_INSPECT_METHOD clotho::recombine::inspection::tag::copy_matching_classify_mismatch
#endif  // USE_INSPECT_ALL

#endif  // USE_VECTOR_SUBSET

#include "clotho/utility/random_generator.hpp"

#ifdef USE_BERNOULLI_RECOMB
#include "clotho/classifiers/bernoulli_classifier_generator.hpp"
#define RECOMBTYPE clotho::classifiers::bernoulli_classifier< rng_type, double, bool>
#else
#include "clotho/classifiers/region_classifier.hpp"
#define RECOMBTYPE clotho::classifiers::region_classifier< allele_type >
#endif  // USE_BERNOULLI_RECOMB


#ifdef USE_BLOCK_SIZE_32
#define BLOCK_UNIT_TYPE unsigned int
#define BLOCK_UNIT_SIZE 32
#else
#define BLOCK_UNIT_TYPE unsigned long
#define BLOCK_UNIT_SIZE 64
#endif  // USE_BLOCK_SIZE_32

#endif  // QTL_ESTIMATOR_CONFIG_HPP_
