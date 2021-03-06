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
#ifndef SPACE_KERNEL_API_HPP_
#define SPACE_KERNEL_API_HPP_

#include <cuda.h>

template < class SpaceType >
__global__ void _delete_space( SpaceType * space );

template < class SpaceType >
__global__ void _resize_space( SpaceType * aspace, unsigned int N );

template < class SpaceType1, class SpaceType2, class OutputSpaceType >
__global__ void _merge_space( SpaceType1 * a, SpaceType2 * b, OutputSpaceType * out );

template < class SpaceType >
__device__ void _update_space( SpaceType * space1, SpaceType * space2 );

#endif  // SPACE_KERNEL_API_HPP_
