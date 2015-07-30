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
#ifndef COMPUTE_CAPABILITY30_HPP_
#define COMPUTE_CAPABILITY30_HPP_

#include "compute_capability/compute_capability_def.hpp"

template < >
struct compute_capability< 3, 0 > {
    static const unsigned int WARP_SIZE = 32;
    static const unsigned int THREADS_PER_BLOCK = 1024;

    static const unsigned int THREAD_BLOCK_2D_X = WARP_SIZE;
    static const unsigned int THREAD_BLOCK_2D_Y = THREADS_PER_BLOCK / WARP_SIZE;

    static const unsigned int WARP_PER_BLOCK = THREADS_PER_BLOCK / WARP_SIZE;

    static const unsigned int MAX_BLOCKS_X = 65535;
    static const unsigned int MAX_BLOCKS_Y = 65535;
};

#endif  // COMPUTE_CAPABILITY30_HPP_
