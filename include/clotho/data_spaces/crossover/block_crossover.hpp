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

namespace clotho {
namespace genetics {

template < class Classifier, class BlockType >
class block_crossover {
public:
    typedef BlockType block_type;
    typedef Classifier classifier_type;

    typedef clotho::utility::debruijn_bit_walker< block_type >  bit_walker_type;
    typedef clotho::utility::BitHelper< block_type >            bit_helper_type;

    block_crossover( const classifier_type & events ) : m_events( events ) {}

/**
 *
 * @param top_strand - block from current top strand
 * @param bottom_strand - block from current bottom strand
 * @param offset - strand relative block offset
 *
 */
    block_type crossover( block_type top_strand, block_type bottom_strand, unsigned int offset ) {
        block_type hets = top_strand ^ bottom_strand;
        block_type mask = bit_helper_type::ALL_UNSET;   // mask state from m_p1 strand
        offset *= bit_helper_type::BITS_PER_BLOCK;      // scale offset to appropriate allele offset

        if( hets ) {
            do {
                unsigned int idx = bit_walker_type::unset_next_index( hets );

                if( m_events( offset + idx ) ) {
                    mask |= clotho::utility::bit_masks[ idx ];
                }
            } while( hets );
        }

        return  ((top_strand & ~mask) | (bottom_strand & mask));
    }

    virtual ~block_crossover() { }

protected:
    classifier_type m_events;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_BLOCK_CROSSOVER_HPP_
