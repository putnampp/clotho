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
#ifndef CLOTHO_SEGMENT_CROSSOVER_TASK_HPP_
#define CLOTHO_SEGMENT_CROSSOVER_TASK_HPP_

#include "clotho/data_spaces/task/task.hpp"
#include "clotho/data_spaces/crossover/block_crossover.hpp"

namespace clotho {
namespace genetics {

template < class Classifier, class BlockType >
class segment_crossover_task : public block_crossover< Classifier, BlockType >, public task {
public:

    typedef block_crossover< Classifier, BlockType > crossover_type;
    typedef BlockType block_type;

    segment_crossover_task( const Classifier & events, block_type * top, block_type * bottom, block_type * o, unsigned int offset, unsigned int len ) :
        crossover_type( events )
        , m_top_strand( top )
        , m_bottom_strand( bottom )
        , m_offspring( o )
        , m_offset( offset )
        , m_length( len )
    {}

    void operator()() {
        const unsigned int end = m_offset + m_length;

        for( unsigned int i = m_offset; i < end; ++i ) {
            block_type t = m_top_strand[ i ];
            block_type b = m_bottom_strand[ i ];

            block_type o = crossover( t, b, i );

            m_offspring[ i ] = o;
        }
    }

    virtual segment_crossover_task() {}

protected:
    block_type      * m_top_strand, * m_bottom_strand, * m_offspring;
    unsigned int    m_offset, m_length;
};

}   // namespace genetics
}   // namespace clotho
#endif  // CLOTHO_SEGMENT_CROSSOVER_TASK_HPP_

