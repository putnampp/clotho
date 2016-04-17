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
#ifndef SIMULATE_ENGINE_BASE_HPP_
#define SIMULATE_ENGINE_BASE_HPP_

#include "qtl_config.hpp"

#include <vector>
#include <utility>

template < class URNG, class AlleleType, class LogType, class TimerType >
class simulate_engine_base {
public:
    typedef URNG                                                        rng_type;
    typedef AlleleType                                                  allele_type;
    typedef LogType                                                     log_type;
    typedef TimerType                                                   timer_type;

    typedef BLOCK_UNIT_TYPE                                             block_type;

    typedef clotho::utility::random_generator< rng_type, allele_type >  allele_generator;

    typedef SUBSETTYPE< allele_type, block_type >   sequence_type;
    typedef typename sequence_type::pointer         sequence_pointer;
    typedef typename sequence_type::powerset_type   allele_set_type;

    typedef RECOMBTYPE                              classifier_type;
    typedef clotho::utility::random_generator< rng_type, classifier_type >  classifier_generator;

//    typedef std::pair< sequence_pointer, sequence_pointer >                  individual_type;
//    typedef std::vector< individual_type >                                   population_type;
//    typedef typename population_type::iterator                               population_iterator;
};

#endif  // SIMULATE_ENGINE_BASE_HPP_
