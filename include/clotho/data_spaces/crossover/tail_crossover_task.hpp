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

#include "clotho/data_spaces/task/task.hpp"
#include "clotho/data_spaces/crossover/tail_crossover.hpp"

namespace clotho {
namespace genetics {

template < class Classifier, class BlockType, class StrandTail >
class tail_crossover_task : public tail_crossover< Classifier, BlockType, StrandTail >, public task {
public:
    typedef tail_crossover< Classifier, BlockType, StrandTail > crossover_type;
    typedef BlockType                                           block_type;

    tail_crossover_task( const Classifier & events, block_type * strand, block_type * offspring, unsigned int offset, unsigned int len ) :
        crossover_type( events )
        , m_strand( strand )
        , m_offspring( offspring )
        , m_offset( offset )
        , m_length( len )
    {
        // should be able to correct upstream logic
        // to prevent 0 length crossover tasks from being created
        // waste of time and memory
#ifdef DEBUGGING
        assert( m_length > 0 );
#endif  // DEBUGGING
        
    }
 
    void operator()() {
        const unsigned int end = m_offset + m_length;
//#ifdef DEBUGGING
//        BOOST_LOG_TRIVIAL(debug) << "Tail crossover: " << m_strand  << " [" << m_offset << ", " << end << "); " << this->m_cfier.event_count();
//#endif  // DEBUGGING

        for( unsigned int i = m_offset; i < end; ++i ) {
            block_type b = m_strand[ i ];

            block_type o = this->crossover( b, i );

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
