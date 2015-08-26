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

template < class SpaceType >
void allele_space_alloc( SpaceType *& space, size_t N = 0 );

template < class SpaceType >
void allele_space_free( SpaceType * space );

template < class SpaceType >
unsigned int resize_space( SpaceType * space, unsigned int N );

template < class SpaceType >
bool check_free_space( SpaceType * space );

template < class SpaceType >
void update_free_count( SpaceType * space );

#endif  // DEVICE_ALLELE_SPACE_API_HPP_
