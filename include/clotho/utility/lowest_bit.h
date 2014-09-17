#ifndef CLOTHO_LOWEST_BIT_H_
#define CLOTHO_LOWEST_BIT_H_

#include <cstdlib>

struct lowest_set_bit_node {
    unsigned char bit_index, bit_shift_next;
    lowest_set_bit_node * next_ptr;

    lowest_set_bit_node() : bit_index(0), bit_shift_next(0), next_ptr(NULL) {}
    lowest_set_bit_node( const lowest_set_bit_node & other ) : bit_index( other.bit_index ), bit_shift_next( other.bit_shift_next ), next_ptr( other.next_ptr) {}
};

class lowest_bit_256 {
public:
    typedef lowest_set_bit_node value_type;
    typedef unsigned char block_type;

    static const unsigned int block_width = 8;
    static const unsigned int max_values = 256;

    lowest_bit_256();

    const value_type & operator[]( block_type val ) const;

    unsigned char low_bit_index( block_type val ) const;
    unsigned char next_shift_offset( block_type val ) const;

    const value_type * begin() const;

    virtual ~lowest_bit_256();
protected:
    void initialize();

    value_type m_vals[max_values];
};

class lowest_bit_64K {
public:
    typedef lowest_set_bit_node value_type;
    typedef unsigned short block_type;

    static const unsigned int block_width = 16;
    static const unsigned int max_values = 65536;

    lowest_bit_64K();

    const value_type & operator[]( block_type val ) const;

    unsigned char low_bit_index( block_type val ) const;
    unsigned char next_shift_offset( block_type val ) const;

    const value_type * begin() const;

    virtual ~lowest_bit_64K();
protected:
    void initialize();

    value_type m_vals[max_values];
};
#endif  // CLOTHO_LOWEST_BIT_H_
