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
#ifndef LIST_TO_STRING_HPP_
#define LIST_TO_STRING_HPP_

#include <sstream>
#include <string>

template < class Iter >
std::string list_to_string( Iter first, Iter last ) {
    if( first == last ) return std::string("");

    std::ostringstream oss;
    oss << *first;
    while( ++first != last ) {
        oss << ":" << *first;
    }

    return oss.str();
}

#endif  // LIST_TO_STRING_HPP_
