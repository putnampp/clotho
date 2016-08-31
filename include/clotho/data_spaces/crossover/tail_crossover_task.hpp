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
#ifndef CLOTHO_TAIL_CROSSOVER_TASK_HPP_
#define CLOTHO_TAIL_CROSSOVER_TASK_HPP_

#include "clotho/data_spaces/crossover/tail_crossover.hpp"

namespace clotho {
namespace genetics {

template < class EventList, class BlockType, class StrandTail >
class tail_crossover_task : public tail_crossover< EventList, BlockType, StrandTail >, public task {
public:
    typedef tail_crossover< EventList, BlockType, StrandTail > crossover_type;

    tail_crossover_task( EventList * events, block_type * strand, block_type * offspring, unsigned int offset, unsigned int len ) :
        crossover_type( events )
        , m_strand( strand )
        , m_offspring( offspring )
        , m_offset( offset )
        , m_length( len )
    {}
 
    void operator()() {
        const unsigned int end = m_offset + m_length;

        for( unsigned int i = m_offset; i < end; ++i ) {
            block_type b = m_strand[ i ];

            block_type o = evaluate( b, i );

            m_offspring[ i ] = o;
        }
    }

    virtual ~tail_crossover_task() {}

protected:
    block_type * m_strand, * m_offspring;

    unsigned int m_offset, m_length;
};


}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_TAIL_CROSSOVER_TASK_HPP_
