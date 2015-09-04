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
#ifndef SPACE_API_HPP_
#define SPACE_API_HPP_

template < class SpaceType >
void create_space( SpaceType *& space, unsigned int N = 0 );

template < class SpaceType >
void resize_space( SpaceType * space, unsigned int N );

template < class SpaceType >
void delete_space( SpaceType * space );

#endif  // SPACE_API_HPP_
