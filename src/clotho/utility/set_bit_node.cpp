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
