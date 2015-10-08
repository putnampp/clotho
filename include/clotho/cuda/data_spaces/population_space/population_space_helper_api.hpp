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
#ifndef POPULATION_SPACE_HELPER_API_HPP_
#define POPULATION_SPACE_HELPER_API_HPP_


template < class IntType, class OrderTag >
void update_free_space2( device_sequence_space< IntType > * seq_space
                       , device_free_space< IntType, OrderTag > * free_space );


template < class IntType, class OrderTag >
void update_free_map( device_free_space< IntType, OrderTag > * fspace );

template < class AlleleSpaceType, class FreeSpaceType >
void update_fixed_alleles( AlleleSpaceType * alleles, FreeSpaceType * free_space );

#endif  // POPULATION_SPACE_HELPER_API_HPP_
