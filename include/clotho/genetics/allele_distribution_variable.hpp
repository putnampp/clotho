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
