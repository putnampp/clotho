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
#ifndef DEVICE_ALLELE_SPACE_API_HPP_
#define DEVICE_ALLELE_SPACE_API_HPP_

#include "clotho/cuda/data_spaces/data_space_api.hpp"

/*
template < class SpaceType >
void create_allele_space( SpaceType *& space, unsigned int N = 0 );

template < class SpaceType >
void delete_allele_space( SpaceType * space );

template < class SpaceType >
unsigned int resize_space( SpaceType * space, unsigned int N );

template < class SpaceType >
void resize_allele_space( SpaceType * space, unsigned int N );


template < class SpaceType >
bool check_free_space( SpaceType * space );

template < class SpaceType >
void update_free_count( SpaceType * space );

template < class SpaceType >
void merge_allele_space( SpaceType * a, SpaceType * b, SpaceType * output );
*/

#include "clotho/cuda/data_spaces/event_space/device_event_space.hpp"

template < class RealType, class IntType, class OrderTag >
void merge_space( device_allele_space< RealType, IntType, OrderTag > * in_space
                , device_event_space< IntType > * evts
                , device_allele_space< RealType, IntType, OrderTag > * out_space );

#endif  // DEVICE_ALLELE_SPACE_API_HPP_
