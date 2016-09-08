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
#ifndef CLOTHO_BLOCK_MUTATE_HPP_
#define CLOTHO_BLOCK_MUTATE_HPP_

#include "clotho/utility/bit_helper.hpp"

namespace clotho {
namespace genetics {

template < class BlockType >
class block_mutate {
public:

    typedef BlockType block_type;

    typedef clotho::utility::BitHelper< block_type >    bit_helper_type;

    block_mutate( unsigned int offset ) :
        m_mask( bit_helper_type::ALL_UNSET )
    {
        m_mask = bit_helper_type::bit_offset( offset );
    }

    template < class OffsetIterator >
    block_mutate( OffsetIterator first, OffsetIterator last ) :
        m_mask( bit_helper_type::ALL_UNSET )
    {
        while( first != last ) {
            block_type b = bit_helper_type::bit_offset( *first++ );

            assert( (m_mask & b) == bit_helper_type::ALL_UNSET );

            m_mask ^= b;
        }
    }

    block_type  mutate( block_type b ) {
        assert( (b & m_mask) == bit_helper_type::ALL_UNSET );

        return (b ^ m_mask);
    }

    virtual ~block_mutate() {}

protected:
    block_type m_mask;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_BLOCK_MUTATE_HPP_

