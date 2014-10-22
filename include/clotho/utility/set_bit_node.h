#ifndef CLOTHO_SET_BIT_NODE_H_
#define CLOTHO_SET_BIT_NODE_H_

namespace clotho {
namespace utility {

struct set_bit_node {
    unsigned char bit_index, bit_shift_next;
    set_bit_node * next_ptr;
};

struct set_bit_node_array {
    static bool init_array( set_bit_node * first, unsigned int max );
};

}   // namespace utility
}   // namespace clotho
#endif  // CLOTHO_SET_BIT_NODE_H_
