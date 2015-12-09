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
#ifndef MERGE_SPACE_HELPER_HPP_
#define MERGE_SPACE_HELPER_HPP_

template < class OrderTag >
struct merge_execution_config {
    static const unsigned int BLOCK_COUNT = 1;
    static const unsigned int THREAD_COUNT = 1;
};

template < class IntType>
struct merge_execution_config < unit_ordered_tag< IntType > > {
    static const unsigned int BLOCK_COUNT = 1;
    static const unsigned int THREAD_COUNT = unit_ordered_tag< IntType >::OBJECTS_PER_UNIT;
};

#endif  // MERGE_SPACE_HELPER_HPP_
