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
#ifndef TOOLKIT_HPP_
#define TOOLKIT_HPP_

#include "clotho/utility/igenerator.hpp"

#include <boost/property_tree/ptree.hpp>

template < class T >
class toolkit {
public:

    std::shared_ptr< T > create_generator( boost::property_tree::ptree & config );
};

#endif  // TOOLKIT_HPP_
