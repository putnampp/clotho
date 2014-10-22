#ifndef CLOTHO_BIT_WALKER_HPP_
#define CLOTHO_BIT_WALKER_HPP_

#include "clotho/utility/set_bit_node.h"

namespace clotho {
namespace utility {

template < unsigned int B >
class bit_block_walker;

template < >
class bit_block_walker< 8 > : public set_bit_node_array {
public:
    typedef unsigned char block_type;
    static const unsigned int max_nodes = 256;  // 2^8
protected:
    static set_bit_node m_nodes[ max_nodes ];

    static bool init() {
        static bool is_init = init_array( m_nodes, max_nodes );
        return _is_init;
    }
};

template < >
class bit_block_walker< 16 > : public set_bit_node_array {
public:
    typedef unsigned short block_type;
    static const unsigned int max_nodes = 256 * 256;    // 2^16

protected:
    static set_bit_node m_nodes[ max_nodes ];
    static bool init() {
        static bool _is_init = init_array(m_nodes, max_nodes);
        return _is_init;
    }
};

template < unsigned int U >
class sub_block_walker;

template <>
class sub_block_walker< 1 > {

};

template <>
class sub_block_walker< 2 > {

};

template <>
class sub_block_walker< 4 > {
public:
    template < class Block, class OP >
    static void apply( 

};

template < class Block, class SubBlock >
class block_walker {
public:
    typedef Block block_type;
    typedef SubBlock sub_block_type;

    typedef sub_block_walker< sizeof(Block) / sizeof(SubBlock) > sub_walker_type;

    template < class OP >
    static void apply( block_type _bits, OP oper ) {
        sub_walker_type::apply( _bits, oper );
    }
};

}   // namespace utility {
}   // namespace clotho {

#endif  // CLOTHO_BIT_WALKER_HPP_
