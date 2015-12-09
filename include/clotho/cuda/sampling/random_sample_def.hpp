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
#ifndef RANDOM_SAMPLE_DEF_HPP_
#define RANDOM_SAMPLE_DEF_HPP_

#include <cuda.h>

template < class StateType, class SourceSpaceType, class IndexSpaceType >
__global__ void random_sample( StateType * states
                               , SourceSpaceType * src
                               , unsigned int N
                               , IndexSpaceType * events );

#endif  // RANDOM_SAMPLE_DEF_HPP_
