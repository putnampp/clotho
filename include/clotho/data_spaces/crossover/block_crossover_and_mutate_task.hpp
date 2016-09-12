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
#ifndef CLOTHO_BLOCK_CROSSOVER_AND_MUTATE_HPP_
#define CLOTHO_BLOCK_CROSSOVER_AND_MUTATE_HPP_

#include "clotho/data_spaces/crossover/block_crossover.hpp"
#include "clotho/data_spaces/mutation/block_mutate.hpp"
#include "clotho/data_spaces/task/task.hpp"

namespace clotho {
namespace genetics {

template < class Classifier, class BlockType >
class block_crossover_and_mutate_task : public block_crossover< Classifier, BlockType >, public block_mutate< BlockType >, public task {
public:
    typedef block_crossover< Classifier, BlockType >    crossover_type;
    typedef block_mutate< BlockType >                   mutate_type;
    typedef BlockType                                   block_type;

    block_crossover_and_mutate_task( const crossover_type & xover, block_type top, block_type bottom, block_type * o, unsigned int index, const mutate_type & mut ) :
        crossover_type( xover )
        , mutate_type( mut )
        , m_top_strand( top )
        , m_bottom_strand( bottom )
        , m_offspring( o )
        , m_index( index )
    {}

    block_crossover_and_mutate_task( const Classifier & events, block_type top, block_type bottom, block_type * o, unsigned int index, const mutate_type & mut ) :
        crossover_type( events )
        , mutate_type ( mut )
        , m_top_strand( top )
        , m_bottom_strand( bottom )
        , m_offspring( o )
        , m_index( index )
    {}

    template < class OffsetIterator >
    block_crossover_and_mutate_task( const Classifier & events, block_type top, block_type bottom, block_type * o, unsigned int index, OffsetIterator first, OffsetIterator last ) :
        crossover_type( events )
        , mutate_type( first, last )
        , m_top_strand( top )
        , m_bottom_strand( bottom )
        , m_offspring(o)
        , m_index( index )
    {}

    void operator()() {
        block_type o = this->crossover( m_top_strand, m_bottom_strand, m_index );

        o = this->mutate( o );

        m_offspring[ m_index ] = o; 
    }

    virtual ~block_crossover_and_mutate_task() {}

protected:
    block_type m_top_strand, m_bottom_strand;
    block_type * m_offspring;

    // block index is necessary for the crossover algorithm
    unsigned int m_index;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_BLOCK_CROSSOVER_AND_MUTATE_HPP_

