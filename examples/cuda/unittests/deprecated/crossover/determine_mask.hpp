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
#ifndef DETERMINE_MASK_HPP_
#define DETERMINE_MASK_HPP_

#include <algorithm>

template < class DataVector, class ElementIterator >
inline unsigned int determine_mask( DataVector & rec_events, ElementIterator first, ElementIterator last ) {
    unsigned int m = 0, bit = 1;
    while( first != last ) {
        typename DataVector::iterator it = std::upper_bound( rec_events.begin(), rec_events.end(), *first );
        m |= (( rec_events.begin() - it ) & 1) * bit;
        bit <<= 1;
        ++first;
    }
    return m;
}
#endif  // DETERMINE_MASK_HPP_
