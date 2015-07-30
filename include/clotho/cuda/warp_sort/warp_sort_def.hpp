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
#ifndef WARP_SORT_DEF_HPP_
#define WARP_SORT_DEF_HPP_

namespace clotho {
namespace cuda {

/**
 * warp_sort: shuffles values in local registers within a warp
 * such that values are in ascending order with the laneid (tid % WARP_WIDTH).
 *
 * Performs a bitonic sort using registers of threads within a warp.
 * The initial version of the algorithm is provided in:
 * 
 * @inproceedings{Demouth2013,
 * author = {Julien Demouth},
 * title = {Shuffle: Tips and Tricks},
 * year = {2013},
 * url = {http://on-demand.gputechconf.com/gtc/2013/presentations/S3174-Kepler-Shuffle-Tips-Tricks.pdf},
 * }
 *
 * @param WW - warp width
 */
template < unsigned int WW >
struct warp_sort;

}   // namespace cuda
}   // namespace clotho

#endif  // WARP_SORT_DEF_HPP_
