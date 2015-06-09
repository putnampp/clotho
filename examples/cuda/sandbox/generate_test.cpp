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
#include <iostream>

#include "generatable_type_trait.hpp"

struct A {
    int operator()() { return 0; }
};

//struct B { };

int main() {
    std::cout << std::boolalpha;

    std::cout << "Is A generatable? " << is_generatable< A >::value << "\n";

      // fails at compile as expected
//    std::cout << "Is B generatable? " << is_generatable< B >::value << "\n";

    return 0;
}
