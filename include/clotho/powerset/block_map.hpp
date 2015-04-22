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
#ifndef BLOCK_MAP_HPP_
#define BLOCK_MAP_HPP_

#include "clotho/powerset/bit_range.hpp"
#include "clotho/powerset/normalized_key.hpp"

#include "clotho/powerset/next_available_order_tag.hpp"

namespace clotho {
namespace powersets {

template < class Element, class Block, class Tag = next_available_order_tag >
struct block_map {
    typedef Element element_type;
    typedef Block   size_type;

    static const unsigned int       bits_per_block = sizeof( size_type ) * 8;

    size_type operator()( const element_type & elem ) {
        return (size_type)-1;
    }
};

template < class Element, class Block >
struct block_map< Element, Block, normalized_key< Element > > : public bit_range< Block, typename normalized_key< Element >::range_type > {
    typedef Element element_type;
    typedef Block   size_type;

    typedef normalized_key< Element > keyer;

    typedef bit_range< Block, typename normalized_key< Element >::range_type > bit_range_type;

//    static const unsigned int       bits_per_block = sizeof( size_type ) * 8;
//    constexpr static const double   width_per_bin = ((keyer::range) / (double) bits_per_block);

    inline size_type operator()( const element_type & elem ) {
        static keyer k;
        return k(elem) * bit_range_type::bits_per_block;
    }
};

}   // namespace powersets
}   // namespace clotho

#endif  // BLOCK_MAP_HPP_
