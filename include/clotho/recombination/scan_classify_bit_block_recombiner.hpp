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
#ifndef SCAN_CLASSIFY_BIT_BLOCK_RECOMBINER_HPP_
#define SCAN_CLASSIFY_BIT_BLOCK_RECOMBINER_HPP_

#include "clotho/recombination/bit_block_recombiner_def.hpp"
#include "clotho/recombination/inspect_methods.hpp"

#include "clotho/utility/bit_masks.hpp"

namespace clotho {
namespace recombine {
namespace walker {
namespace tag {

/**
 * Scan and Classify walks the set bits in a block and stores their
 * index (explodes the block into list of bit indices). The list of
 * indices is iterated and the elements at these offsets are classified.
 *
 * This results in:
 *      - a single bit walking loop per 32-bit block
 *      - a simple loop over offsets
 *      - switch over the bit indices
 *
 * Thought process:
 *      - Bit scanning could more streamline (less context switching)
 *      - Requires more memory (64 * 4 = 256 Bytes)
 *      - Simpler classification loop (indices will be ordered as a result of scan)
 *      - Switch logic allows mask bits sets to be defined inline rather than computed
 *
 *  Concern:
 *      - Doubles the number of iterative steps
 *      - Bit mask computation dependent upon index (adds bit shifting step)
 */
struct scan_and_classify {};

}   // namespace tag
}   // namespace walker
}   // namespace recombine
}   // namespace clotho

namespace clotho {
namespace recombine {

template < class Classifier, class InspectMethodTag >
class bit_block_recombiner< Classifier, InspectMethodTag, clotho::recombine::walker::tag::scan_and_classify > {
public:
    typedef Classifier classifier_type;

    bit_block_recombiner( const classifier_type & cfier ) : m_cfier( cfier ) {}

    template < class Block, class ElementIterator >
    Block operator()( const Block b0, const Block b1, const ElementIterator first ) {
        unsigned int count = clotho::utility::scan_set_bits( InspectMethodTag::select(b0, b1), indices );
        unsigned int * idx = indices;

        Block res = (Block)0;
        while( count-- ) {
            unsigned int tmp = *idx++;
            if( m_cfier( *(first + tmp) ) ) {
                res |= ((Block)1 << tmp);
            }
        }

        Block rec = ((b0 & res) | (b1 & ~res));
        return rec;
    }

protected:
    classifier_type m_cfier;
    unsigned int indices[ 64 ];
};

}   // namespace recombine
}   // namespace clotho
#endif  // SCAN_CLASSIFY_BIT_BLOCK_RECOMBINER_HPP_
