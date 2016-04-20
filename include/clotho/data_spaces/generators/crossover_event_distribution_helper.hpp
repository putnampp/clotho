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
#ifndef CLOTHO_CROSSOVER_EVENT_DISTRIBUTION_HELPER_HPP_
#define CLOTHO_CROSSOVER_EVENT_DISTRIBUTION_HELPER_HPP_

#include <boost/random/poisson_distribution.hpp>

namespace clotho {
namespace genetics {

template < class T >
struct crossover_event_distribution_helper;

template < >
struct crossover_event_distribution_helper< double > {
    typedef unsigned int    IntType;

    typedef boost::random::poisson_distribution< IntType, double > type;
};

template < >
struct crossover_event_distribution_helper< float > {
    typedef unsigned int    IntType;
    typedef boost::random::poisson_distribution< IntType, float > type;
};


}   // namespace genetics
}   // namespace cloth

#endif  // CLOTHO_CROSSOVER_EVENT_DISTRIBUTION_HELPER_HPP_
