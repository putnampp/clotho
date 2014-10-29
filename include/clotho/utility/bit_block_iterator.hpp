#ifndef BIT_BLOCK_ITERATOR_HPP_
#define BIT_BLOCK_ITERATOR_HPP_

#include <type_traits>

namespace clotho {
namespace utility {

struct linear_iterator_tag {};

template < class Block, class Tag = linear_iterator_tag, class Enabled = void >
class bit_block_iterator;

template < class Block >
class bit_block_iterator < Block, linear_iterator_tag, typename std::enable_if< std::is_integral< Block >::value >::type > {
public:
    typedef Block block_type;

    bit_block_iterator( block_type b = (block_type)0) :
        m_val(b)
        , m_index(0)
    {
        if( m_val & (block_type)1 == 0) 
            next();
    }

    bit_block_iterator( const bit_block_iterator< Block > & rhs ) :
        m_val( rhs.m_val )
        , m_index( rhs.m_index )
    {}

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

    bool operator==( const bit_block_iterator< Block > & rhs ) const {
        return (this->m_val == rhs.m_val);
    }

    bool operator!=(const bit_block_iterator< Block > & rhs ) const {
        return (this->m_val != rhs.m_val);
    }

    virtual ~bit_block_iterator() {}

protected:
    inline void next() {
        do {
            m_val >>= 1;
            ++m_index;
        } while( (m_val & (block_type)1 == 0) );
    }

    block_type m_val;
    unsigned int m_index;
};

}   // namespace utility {
}   // namespace clotho {

#endif  // BIT_BLOCK_ITERATOR_HPP_
