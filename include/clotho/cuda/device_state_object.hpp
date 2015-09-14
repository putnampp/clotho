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
#ifndef DEVICE_GET_STATE_HPP_
#define DEVICE_GET_STATE_HPP_

#include "clotho/utility/state_object.hpp"

template < class ObjectType >
void get_device_object_state( boost::property_tree::ptree & state, ObjectType * obj ) {
    ObjectType local_obj;

    assert( cudaMemcpy( &local_obj, obj, sizeof( ObjectType ), cudaMemcpyDeviceToHost ) == cudaSuccess );

    clotho::utility::get_state( state, local_obj );
}

#endif  // DEVICE_GET_STATE_HPP_
