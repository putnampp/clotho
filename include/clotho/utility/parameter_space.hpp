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
#ifndef PARAMETER_SPACE_HPP_
#define PARAMETER_SPACE_HPP_

#include <boost/property_tree/ptree.hpp>

namespace clotho {
namespace utility {

template < class Configable >
struct parameter_space {
    static void build_parameters( boost::property_tree::ptree & params ) {}
};

}   // namespace utility
}   // namespace clotho

#endif  // PARAMETER_SPACE_HPP_
