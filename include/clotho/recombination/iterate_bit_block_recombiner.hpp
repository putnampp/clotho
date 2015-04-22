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
#ifndef ITERATE_BIT_BLOCK_RECOMBINER_HPP_
#define ITERATE_BIT_BLOCK_RECOMBINER_HPP_

#include "clotho/recombination/bit_block_recombiner_def.hpp"
#include "clotho/recombination/inspect_methods.hpp"

namespace clotho {
namespace recombine {
namespace walker {
namespace tag {

struct iterate_and_classify {};

}   // namespace tag
}   // namespace walker
}   // namespace recombine
}   // namespace clotho

namespace clotho {
namespace recombine {

template < class Classifier, class InspectMethodTag >
class bit_block_recombiner< Classifier, InspectMethodTag, clotho::recombine::walker::tag::iterate_and_classify > {
public:
    typedef Classifier classifier_type;

    bit_block_recombiner( const classifier_type & cfier ) : m_cfier( cfier ) {}

    template < class Block, class ElementIterator >
    Block operator()( const Block b0, const Block b1, const ElementIterator first ) {
        typedef clotho::utility::bit_block_iterator< Block, clotho::utility::tag::linear_iterator_tag > iterator;

        Block mask = (Block)0;
        iterator bit_it( InspectMethodTag::select(b0, b1) ), bit_end;
        while( bit_it != bit_end ) {
            unsigned int idx = (*bit_it++);
            if( m_cfier( *(first + idx) ) ) {
                mask |= ((Block)1 << idx );
            }
        }

        Block res = ((b0 & mask) | (b1 & ~ mask) );
        return res;

    }

protected:
    classifier_type m_cfier;
};

}   // namespace recombine
}   // namespace clotho
#endif  // ITERATE_BIT_BLOCK_RECOMBINER_HPP_
