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
#include "clotho/utility/set_bit_node.h"

#include <iostream>

namespace clotho {
namespace utility {


bool node_array_initializer::init_array( set_bit_node * first, unsigned int max ) {
    unsigned int val = max - 1;
    while( val ) {
        unsigned int idx = 0, tmp = val;

        while( tmp ) {
            if( tmp & 1 ) {
                first[val].bit_index = idx++;
                tmp >>= 1;
                break;
            }
            tmp >>= 1;
            ++idx;
        }

        if( tmp ) {
            do {
                if( tmp & 1 ) {
                    first[ val ].bit_shift_next = idx;
                    first[ val ].next = tmp;
                    break;
                }
                tmp >>= 1;
                ++idx;
            } while( tmp );
        } else {
            first[val].bit_shift_next = idx;
            first[val].next = 0;
        }
        --val;
    }

    return true;
}

}   // namespace utility
}   // namespace clotho
