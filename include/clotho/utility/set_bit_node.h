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
#ifndef CLOTHO_SET_BIT_NODE_H_
#define CLOTHO_SET_BIT_NODE_H_

namespace clotho {
namespace utility {

struct set_bit_node {
    unsigned char bit_index, bit_shift_next;
    unsigned int next;
    set_bit_node( unsigned char i = 0, unsigned char s = 0, unsigned char n = 0 ) :
        bit_index(i), bit_shift_next(s), next(n) {}
};

struct node_array_initializer {
    static bool init_array( set_bit_node * first, unsigned int max );
};

}   // namespace utility
}   // namespace clotho

#endif  // CLOTHO_SET_BIT_NODE_H_
