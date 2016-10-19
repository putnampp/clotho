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
#ifndef DEBRUIJN_BIT_BLOCK_ITERATOR2_HPP_
#define DEBRUIJN_BIT_BLOCK_ITERATOR2_HPP_

#include "clotho/utility/bit_block_iterator_def.hpp"
#include "clotho/utility/bit_masks.hpp"

namespace clotho {
namespace utility {
namespace tag {

struct debruijn_iterator_tag2 {};

}   // namespace tag
}   // namespace utility
}   // namespace clotho

namespace clotho {
namespace utility {

template < class Block >
class bit_block_iterator < Block, clotho::utility::tag::debruijn_iterator_tag2, typename std::enable_if< std::is_integral< Block >::value >::type > {
public:
    typedef bit_block_iterator< Block, clotho::utility::tag::debruijn_iterator_tag2, void > self_type;
    typedef Block block_type;

    static const unsigned int bits_per_block = sizeof( block_type ) * 8;

    bit_block_iterator( block_type b = (block_type)0) :
        m_val(b),
        m_lsb(0) 
    { }

    bit_block_iterator( const self_type & rhs ) :
        m_val( rhs.m_val )
        , m_lsb( rhs.m_lsb ) {
    }

    bool hasNext() const {
        return m_val;
    }

    unsigned int next() {
        unsigned int lo = m_val;
        if( !lo ) {
            m_lsb += 32;
            m_val >>= 32;
            lo = m_val;
        }

        lo = LEAST_SIG_BIT( lo );
        lo = DEBRUIJNBIT_HASH_LOOKUP( lo );

        m_val >>= lo;
        m_val ^= (block_type) 1;

        m_lsb += lo;

        return m_lsb;
    }

    bool operator==( const self_type & rhs ) const {
        return (this->m_val == rhs.m_val);
    }

    bool operator!=(const self_type & rhs ) const {
        return (this->m_val != rhs.m_val );
    }

    void reset( block_type b ) {
        m_val = b;
        m_lsb = 0;
    }

    virtual ~bit_block_iterator() {}

protected:
    block_type m_val, m_lsb;
};

}   // namespace utility
}   // namespace clotho
#endif  // DEBRUIJN_BIT_BLOCK_ITERATOR2_HPP_
