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
#ifndef VALDIATE_ORDER_HPP_
#define VALDIATE_ORDER_HPP_

template < class Iter1, class Iter2 >
bool validate_order( Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2 ) {
    const double epsilon = 0.000001;
    bool eq = true;
    while( eq ) {
        if( first1 == last1 ) {
            eq = (first2 == last2);
            break;
        } else if( first2 == last2 ) {
            eq = false;
            break;
        }

        eq = ( abs(*first1 - *first2) < epsilon );
        ++first1; ++first2;
    }
    return eq;
}

#endif  // VALDIATE_ORDER_HPP_
