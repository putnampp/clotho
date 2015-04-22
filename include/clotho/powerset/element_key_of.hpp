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
#ifndef ELEMENT_KEY_OF_HPP_
#define ELEMENT_KEY_OF_HPP_


namespace clotho {
namespace powersets {

template < class E >
struct element_key_of {
    typedef E key_type;

    key_type operator()( const E & e ) {
        return e;
    }

    static key_type get_key( const E & e ) {
        return e;
    }
};

}   // namespace powersets
}   // namespace cl

#endif  // ELEMENT_KEY_OF_HPP_
