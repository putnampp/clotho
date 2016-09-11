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
#ifndef CLOTHO_BLOCK_MUTATE_TASK_HPP_
#define CLOTHO_BLOCK_MUTATE_TASK_HPP_

#include "clotho/data_spaces/task/task.hpp"
#include "clotho/data_spaces/mutation/block_mutate.hpp"

namespace clotho {
namespace genetics {

template < class BlockType >
class block_mutate_task : public block_mutate< BlockType >, public task {
public:
    typedef block_mutate< BlockType >   mutate_type;
    typedef BlockType                   block_type;

    block_mutate_task( block_type t, block_type * res, const mutate_type & mut ) :
        mutate_type( mut )
    {}

    void operator()() {
        *m_result = mutate( m_orig );
    }

    virtual ~block_mutate_task() {}
protected:
    block_type m_orig;
    block_type * m_result;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_BLOCK_MUTATE_TASK_HPP_

