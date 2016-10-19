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
#ifndef CLOTHO_BLOCK_CROSSOVER_HPP_
#define CLOTHO_BLOCK_CROSSOVER_HPP_

#include "clotho/utility/bit_helper.hpp"
#include "clotho/utility/debruijn_bit_walker.hpp"
//#include "clotho/utility/bit_block_iterator.hpp"

namespace clotho {
namespace genetics {

template < class Classifier, class BlockType >
class block_crossover {
public:
    typedef BlockType block_type;
    typedef Classifier classifier_type;

    typedef clotho::utility::debruijn_bit_walker< block_type >  bit_walker_type;
    typedef clotho::utility::BitHelper< block_type >            bit_helper_type;

   // typedef clotho::utility::bit_block_iterator< block_type, clotho::utility::tag::debruijn_iterator_tag2 > bit_iterator;

    block_crossover( const classifier_type & events ) :
        m_cfier( events ) 
    { }

/**
 *
 * @param top_strand - block from current top strand
 * @param bottom_strand - block from current bottom strand
 * @param offset - strand relative block offset
 *
 */
    inline block_type crossover( const block_type & top_strand, const block_type & bottom_strand, const unsigned int OFFSET ) {
        block_type hets = top_strand ^ bottom_strand;
        block_type mask = bit_helper_type::ALL_UNSET;   // mask state from m_p1 strand
//        offset *= bit_helper_type::BITS_PER_BLOCK;      // scale offset to appropriate allele offset

        unsigned int idx = 0;
        while( hets ) {
            idx += bit_walker_type::next_and_shift( hets );

            if( m_cfier( OFFSET + idx ) ) {
                // the state of the corresponding allele is defined by the bottom strand
//                mask |= clotho::utility::bit_masks[ idx ];
                mask |= (((block_type)1) << idx);
            }
        }

        return  ((top_strand & ~mask) | (bottom_strand & mask));
    }

    virtual ~block_crossover() { }

protected:
    classifier_type m_cfier;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_BLOCK_CROSSOVER_HPP_
