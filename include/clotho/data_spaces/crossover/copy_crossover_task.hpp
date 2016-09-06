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
#ifndef CLOTHO_COPY_CROSSOVER_TASK_HPP_
#define CLOTHO_COPY_CROSSOVER_TASK_HPP_

#include "clotho/data_spaces/task/task.hpp"

namespace clotho {
namespace genetics {

template < class BlockType >
class copy_crossover_task : public task {
public:
    typedef copy_crossover_task< BlockType > self_type;

    typedef BlockType block_type;

    copy_crossover_task( block_type * strand, block_type * offspring, unsigned int len ) :
        m_strand( strand )
        , m_offspring( offspring )
        , m_length( len )
    {}

    void operator()() {
        memcpy( m_offspring, m_strand, m_length * sizeof( block_type ) );
    }

    virtual ~copy_crossover_task() {}

protected:
    block_type  * m_strand, * m_offspring;
    unsigned int m_length;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_COPY_CROSSOVER_TASK_HPP_
