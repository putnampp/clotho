#ifndef BIT_BLOCK_ITERATOR_HPP_
#define BIT_BLOCK_ITERATOR_HPP_

#include <iostream>
#include <type_traits>
#include "clotho/utility/bit_walker.hpp"

namespace clotho {
namespace utility {

struct linear_iterator_tag {};
struct walk_iterator_tag {};

template < class Block, class Tag = linear_iterator_tag, class Enabled = void >
class bit_block_iterator;

template < class Block >
class bit_block_iterator < Block, linear_iterator_tag, typename std::enable_if< std::is_integral< Block >::value >::type > {
public:
    typedef bit_block_iterator< Block, linear_iterator_tag, void > self_type;
    typedef Block block_type;

    bit_block_iterator( block_type b = (block_type)0) :
        m_val(b)
        , m_index(0) {
        if( (m_val & (block_type)1) == 0)
            next();
    }

    bit_block_iterator( const self_type & rhs ) :
        m_val( rhs.m_val )
        , m_index( rhs.m_index ) {
    }

    bit_block_iterator & operator++() {
        if( m_val )
            next();
        return *this;
    }

    bit_block_iterator operator++(int) {
        bit_block_iterator res(*this);
        this->operator++();
        return res;
    }

    unsigned int operator*() const {
        return m_index;
    }

    bool operator==( const self_type & rhs ) const {
        return (this->m_val == rhs.m_val);
    }

    bool operator!=(const self_type & rhs ) const {
        return (this->m_val != rhs.m_val);
    }

    virtual ~bit_block_iterator() {}

protected:
    inline void next() {
        do {
            m_val >>= 1;
            ++m_index;
        } while( (m_val != 0) && (m_val & (block_type)1) == 0 );
    }

    block_type m_val;
    unsigned int m_index;
};

template < class Block >
class bit_block_iterator < Block, walk_iterator_tag, typename std::enable_if< std::is_integral< Block >::value >::type > {
public:
    typedef bit_block_iterator< Block, walk_iterator_tag, void > self_type;
    typedef Block block_type;

    bit_block_iterator( block_type b = (block_type)0) :
        m_val(b)
        , m_index(0)
        , m_base_index(0)
        , m_shift_next(0) {
        next();
    }

    bit_block_iterator( const self_type & rhs ) :
        m_val( rhs.m_val )
        , m_index( rhs.m_index )
        , m_base_index( rhs.m_base_index )
        , m_shift_next(rhs.m_shift_next) {
    }

    bit_block_iterator & operator++() {
        m_base_index += m_shift_next;
        m_val >>= m_shift_next;
        next();
        return *this;
    }

    bit_block_iterator operator++(int) {
        bit_block_iterator res(*this);
        this->operator++();
        return res;
    }

    unsigned int operator*() const {
        return m_index;
    }

    bool operator==( const self_type & rhs ) const {
        return (this->m_val == rhs.m_val);
    }

    bool operator!=(const self_type & rhs ) const {
        return (this->m_val != rhs.m_val);
    }

    void reset( block_type b ) {
        m_val = b;
        m_index = 0;
        m_base_index = 0;
        m_shift_next = 0;

        next();
    }

    virtual ~bit_block_iterator() {}

protected:
    typedef unsigned short sub_block_type;
    typedef clotho::utility::bit_block_walker< sizeof( sub_block_type ) * 8 > walker_type;

    inline void next() {
        if( m_val ) {
            sub_block_type v = (sub_block_type) m_val;
            while(!v) {
                m_val >>= walker_type::bits_per_block;
                m_base_index += walker_type::bits_per_block;
                v = (sub_block_type) m_val;
            }
            set_bit_node * n = walker_type::getNode() + v;

            m_index = m_base_index + n->bit_index;
            m_shift_next =  n->bit_shift_next;
        }
    }

    block_type m_val;
    unsigned int m_index, m_base_index, m_shift_next;
};
}   // namespace utility {
}   // namespace clotho {

#endif  // BIT_BLOCK_ITERATOR_HPP_
