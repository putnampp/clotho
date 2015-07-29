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
#ifndef COMPUTE_CAPABILITIES_HPP_
#define COMPUTE_CAPABILITIES_HPP_

template < unsigned int V >
struct compute_capability;

template < >
struct compute_capability< 3 > {
    static const unsigned int WARP_SIZE = 32;
    static const unsigned int THREADS_PER_BLOCK = 1024;

    static const unsigned int THREAD_BLOCK_2D_X = WARP_SIZE;
    static const unsigned int THREAD_BLOCK_2D_Y = THREADS_PER_BLOCK / WARP_SIZE;

    static const unsigned int WARP_PER_BLOCK = THREADS_PER_BLOCK / WARP_SIZE;
};

#endif  // COMPUTE_CAPABILITIES_HPP_
