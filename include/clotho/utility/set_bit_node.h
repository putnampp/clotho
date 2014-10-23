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
