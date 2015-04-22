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
#ifndef TRAIT_HELPER_HPP_
#define TRAIT_HELPER_HPP_

template < class Sequence >
struct trait_helper {
    typedef typename Sequence::value_type                        element_type;
    typedef typename element_type::weight_type    weight_type;
};

template < class Sequence >
struct trait_helper< std::shared_ptr< Sequence > > : public trait_helper< Sequence > {
    typedef typename trait_helper< Sequence >::element_type element_type;
    typedef typename trait_helper< Sequence >::trait_weight_type weight_type;
};

#include "clotho/powerset/variable_subset.hpp"

template < class E, class B, class BM, class EK >
struct trait_helper< clotho::powersets::variable_subset< E, B, BM, EK > > {
    typedef E                               element_type;
    typedef typename E::weight_type         weight_type;
};

#endif  // TRAIT_HELPER_HPP_
