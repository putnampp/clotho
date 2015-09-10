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
#ifndef GENERATE_MUTATION_KERNEL_API_HPP_
#define GENERATE_MUTATION_KERNEL_API_HPP_

template < class StateType, class FreeSpaceType, class EventSpaceType, class AlleleSpaceType >
__global__ void _generate_mutation_kernel( StateType * states, FreeSpaceType * fspace, EventSpaceType * events, AlleleSpaceType * alleles );

#endif  // GENERATE_MUTATION_KERNEL_API_HPP_
