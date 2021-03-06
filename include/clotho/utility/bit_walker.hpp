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
#ifndef CLOTHO_BIT_WALKER_HPP_
#define CLOTHO_BIT_WALKER_HPP_

#include <iostream>
#include <cassert>
#include <type_traits>

#include "clotho/utility/set_bit_node.h"

namespace clotho {
namespace utility {

template < unsigned int B, class E = void >
class bit_block_walker;

template < unsigned int B >
class bit_block_walker< B, typename std::enable_if< B <= 64 >::type > {
public:
    typedef size_t              walkable_block_type;
    static const unsigned int   bits_per_block = B;

    template < class OP >
    void operator()( walkable_block_type _bits, OP & oper, unsigned int offset = 0 ) {
        if( !_bits ) return;

        unsigned int i = B, j = offset;
        while( i-- ) {
            if( _bits & 1 ) {
                oper( j );
            }
            _bits >>= 1;
            ++j;
        }
    }
};

template < >
class bit_block_walker< 8, void > {
public:
    typedef unsigned char walkable_block_type;
    static const unsigned int bits_per_block = 8;
    static const unsigned int max_nodes = 256;  // 2^8

    static set_bit_node (&getNode())[ max_nodes ] {
        static bool i = node_array_initializer::init_array( m_nodes, max_nodes );
        assert( i );
        return m_nodes;
    }

    template < class OP >
    void operator()( walkable_block_type _bits, OP & oper, unsigned int offset = 0 ) {
        if( !_bits ) return;

        set_bit_node * v = getNode() + _bits;
        do {
            oper( offset + v->bit_index );
            if( v->next == 0 ) break;

            offset += v->bit_shift_next;
            v = getNode() + v->next;
        } while( true );
    }
protected:
    static set_bit_node m_nodes[ max_nodes ];
};

set_bit_node bit_block_walker< 8 >::m_nodes[ bit_block_walker< 8 >::max_nodes ];

template < >
class bit_block_walker< 16, void > {
public:
    typedef unsigned short walkable_block_type;
    static const unsigned int bits_per_block = 16;
    static const unsigned int max_nodes = 256 * 256;    // 2^16

    static set_bit_node (&getNode())[ max_nodes ] {
        static bool i = node_array_initializer::init_array(m_nodes, max_nodes);
//        assert( i );
        return m_nodes;
    }

    template < class OP >
    void operator()( walkable_block_type _bits, OP & oper, unsigned int offset = 0 ) {
        if( !_bits ) return;

        set_bit_node * v = getNode() + _bits;
        do {
            oper( offset + v->bit_index );
            if( v->next == 0 ) break;

            offset += v->bit_shift_next;
            v = getNode() + v->next;
        } while( true );
    }

protected:
    static set_bit_node m_nodes[ max_nodes ];
};

set_bit_node bit_block_walker< 16 >::m_nodes[ bit_block_walker< 16 >::max_nodes ];

template < class Block, class SubBlock, class BitWalker = bit_block_walker< sizeof( SubBlock) * 8> >
class block_walker {
public:
    typedef Block block_type;
    typedef SubBlock sub_block_type;

    static const unsigned int bits_per_block = sizeof( block_type ) * 8;
    static const unsigned int bits_per_subblock = sizeof( sub_block_type ) * 8;

    typedef BitWalker bit_walker_type;

    static const unsigned int subblocks_per_block = (bits_per_block / bit_walker_type::bits_per_block);

    template < class OP >
    static void apply( block_type _bits, OP & oper ) {
        if(! _bits ) return;

        unsigned int i = subblocks_per_block, offset = 0;
        bit_walker_type bw;
        while( i-- ) {
            typename bit_walker_type::walkable_block_type sb = (typename bit_walker_type::walkable_block_type)_bits;
            bw( sb, oper, offset );
            _bits >>= bits_per_subblock;

            if( ! _bits ) break;

            offset += bits_per_subblock;
        }
    }
};

}   // namespace utility {
}   // namespace clotho {

#endif  // CLOTHO_BIT_WALKER_HPP_
