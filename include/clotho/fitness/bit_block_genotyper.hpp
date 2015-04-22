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
#ifndef CLOTHO_BIT_BLOCK_GENOTYPER_HPP_
#define CLOTHO_BIT_BLOCK_GENOTYPER_HPP_

#include <type_traits>

namespace clotho {
namespace fitness {

template < class Block, class Enabled = void >
class bit_block_genotyper;

template < class Block >
class bit_block_genotyper< Block, typename std::enable_if< std::is_integral< Block >::value >::type >  {
public:
    typedef Block block_type;

    static block_type   get_all_heterozygous( const block_type b0, const block_type b1 ) {
        return (b0 ^ b1);
    }

    static block_type get_all_homozygous( const block_type b0, const block_type b1 ) {
        return (~(b0 ^ b1));
    }

    static block_type   get_alt_homozygous( const block_type b0, const block_type b1 ) {
        return (b0 & b1);
    }

    static block_type   get_ref_homozygous( const block_type b0, const block_type b1 ) {
        return (~(b0 | b1));
    }
};

}   // namespace fitness {
}   // namespace clotho {

#endif  // CLOTHO_BIT_BLOCK_GENOTYPER_HPP_
