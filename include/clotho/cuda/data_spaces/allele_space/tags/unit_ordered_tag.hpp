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
#ifndef UNIT_ORDERED_TAG_HPP_
#define UNIT_ORDERED_TAG_HPP_


template < class UnitType = unsigned int >
struct unit_ordered_tag {
    typedef UnitType    unit_type;
    static const unsigned int OBJECTS_PER_UNIT = (sizeof( unit_type ) * 8);
};

#endif  // UNIT_ORDERED_TAG_HPP_
