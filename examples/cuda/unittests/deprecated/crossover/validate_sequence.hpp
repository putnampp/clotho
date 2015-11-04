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
#ifndef SEQUENCE_VALIDATE_HPP_
#define SEQUENCE_VALIDATE_HPP_

template < class Value >
inline bool is_empty_value( const Value & v ) {
    return (v == 0);
}

template < class Iter >
inline bool validate_empty_sequence( Iter first, Iter last ) {
    bool _empty = true;
    while( _empty && first != last ) {
        _empty = is_empty_value( *first );
        ++first;
    }
    return _empty;
}

#endif  // SEQUENCE_VALIDATE_HPP_
