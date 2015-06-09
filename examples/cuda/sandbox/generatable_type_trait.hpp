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
#ifndef GENERATABLE_TRAIT_TYPE_HPP_
#define GENERATABLE_TRAIT_TYPE_HPP_

#include <type_traits>

template < class T >
struct is_generatable : std::is_member_function_pointer< decltype(&T::operator()) > {};

#endif  // GENERATABLE_TRAIT_TYPE_HPP_
