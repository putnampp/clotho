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
#ifndef STATE_OBJECT_HPP_
#define STATE_OBJECT_HPP_

#include <boost/property_tree/ptree.hpp>

namespace clotho {
namespace utility {

struct iStateObject {
    virtual void get_state( boost::property_tree::ptree & state ) = 0;
};

//template < class Object >
//void get_state( boost::property_tree::ptree & s, const Object & obj );

template < class ObjectType >
struct state_getter {
    void operator()( boost::property_tree::ptree & s, const ObjectType & obj );
};

}   // namespace utility
}   // namespace clotho

#endif  // STATE_OBJECT_HPP_
