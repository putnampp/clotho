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
#ifndef CLOTHO_TAIL_CROSSOVER_HPP_
#define CLOTHO_TAIL_CROSSOVER_HPP_

#include "clotho/data_spaces/crossover/block_crossover.hpp"

namespace clotho {
namespace genetics {

template < class Classifier, class BlockType, class StrandType >
class tail_crossover;

struct top_strand_tail {};
struct bottom_strand_tail {};

template < class Classifier, class BlockType >
class tail_crossover< Classifier, BlockType, top_strand_tail > : public block_crossover< Classifier, BlockType > {
public:
    typedef block_crossover< Classifier, BlockType > crossover_type;

    tail_crossover( const Classifier & events ) : crossover_type( events ) {}

    void crossover( block_type b, unsigned int offset ) {
        crossover( b, crossover_type::bit_helper_type::ALL_UNSET, offset );
    }

    virtual ~tail_crossover() {}
};


template < class Classifier, class BlockType >
class tail_crossover< Classifier, BlockType, bottom_strand_tail > : public block_crossover< Classifier, BlockType > {
public:
    typedef block_crossover< Classifier, BlockType > crossover_type;

    tail_crossover( const Classifier & events ) : crossover_type( events ) {}

    void crossover( block_type b, unsigned int offset ) {
        crossover(crossover_type::bit_helper_type::ALL_UNSET, b, offset );
    }

    virtual ~tail_crossover() {}
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_TAIL_CROSSOVER_HPP_
