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
#ifndef CLOTHO_TASK_GROUP_HPP_
#define CLOTHO_TASK_GROUP_HPP_

#include "clotho/data_spaces/task/task.hpp"

#include <memory>

namespace clotho {
namespace genetics {

class TaskGroup {
public:

    void add( std::shared_ptr< task >  t );

    virtual ~TaskGroup() {}
protected:
    vector< std::shared_ptr< task > > m_tasks;
};

}   // namespace genetics
}   // namespace clotho
#endif  // CLOTHO_TASK_GROUP_HPP_
