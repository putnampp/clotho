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
#ifndef CLOTHO_POSITION_DISTRIBUTION_HELPER_HPP_
#define CLOTHO_POSITION_DISTRIBUTION_HELPER_HPP_

#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int_distribution.hpp>

namespace clotho {
namespace genetics {

template < class T >
struct position_distribution_helper;


template < >
struct position_distribution_helper< double > {
    typedef boost::random::uniform_01< double > type;

    static constexpr double MIN = 0.0;
    static constexpr double MAX = 1.0;
};

template < >
struct position_distribution_helper< float > {
    typedef boost::random::uniform_01< float > type;

    static constexpr float MIN = 0.0;
    static constexpr float MAX = 1.0;
};

template < >
struct position_distribution_helper< unsigned int > {
    typedef boost::random::uniform_int_distribution< unsigned int > type;

    static const unsigned int MIN = 0;
    static const unsigned int MAX = 0xFFFFFFFF;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_POSITION_DISTRIBUTION_HELPER_HPP_
