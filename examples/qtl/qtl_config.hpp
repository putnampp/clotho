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
#ifndef QTL_ESTIMATOR_CONFIG_HPP_
#define QTL_ESTIMATOR_CONFIG_HPP_

#include <boost/mpl/aux_/na.hpp>

typedef double qtl_real_type;

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

#include "clotho/classifiers/binary_classifier.hpp"

#ifdef USE_BERNOULLI_RECOMB
#include "clotho/classifiers/bernoulli_classifier_generator.hpp"
#include "clotho/classifiers/tags/is_nonzero_tag.hpp"
#define RECOMBTYPE  clotho::classifiers::binary_classifier< clotho::classifiers::bernoulli_classifier< rng_type, qtl_real_type >, clotho::classifiers::tags::is_nonzero_tag >
#else
#include "clotho/classifiers/region_classifier.hpp"
#include "clotho/classifiers/tags/is_even_tag.hpp"
#define RECOMBTYPE  clotho::classifiers::binary_classifier< clotho::classifiers::region_classifier< allele_type >, clotho::classifiers::tags::is_even_tag >
#endif  // USE_BERNOULLI_RECOMB

#ifdef USE_MUTATE_AND_RECOMBINE
#define REPRODUCTION_METHOD_TAG mutate_recombine_tag
#else
#define REPRODUCTION_METHOD_TAG recombine_mutate_tag
#endif  //MUTATE_AND_RECOMBINE

#ifdef USE_ASSORTATIVE_SELECTOR
#define IND_SELECT assortative_selector
#else
#define IND_SELECT individual_selector
#endif

#ifdef USE_BLOCK_SIZE_32
#define BLOCK_UNIT_TYPE unsigned int
#define BLOCK_UNIT_SIZE 32
#else
#define BLOCK_UNIT_TYPE unsigned long
#define BLOCK_UNIT_SIZE 64
#endif  // USE_BLOCK_SIZE_32

#endif  // QTL_ESTIMATOR_CONFIG_HPP_
