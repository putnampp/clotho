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
#ifndef CLOTHO_BIT_RANGE_HPP_
#define CLOTHO_BIT_RANGE_HPP_

#include "clotho/powerset/key_range.hpp"

namespace clotho {
namespace powersets {

template < class Block, class Range >
struct bit_range : public key_range < Range > {
    typedef Block block_type;

    static const unsigned int       bits_per_block = sizeof( block_type ) * 8;
    constexpr static const double   width_per_bin = ((key_range< Range >::range) / (double) bits_per_block);
};

}   // namespace powersets {
}   // namespace clotho {

#endif  // CLOTHO_BIT_RANGE_HPP_
