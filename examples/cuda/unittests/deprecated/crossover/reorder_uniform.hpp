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
#ifndef REORDER_UNIFORM_HPP_
#define REORDER_UNIFORM_HPP_

#include <cmath>

struct ascending_order {
    template < class RealType >
    inline RealType operator()( RealType v ) {
        return (1.0 - v);
    }
};

struct descending_order {
    template < class RealType >
    inline RealType operator()( RealType v ) {
        return v;
    }
};

/**
 * Do not like this method as it performs two natural log calculations
 */
template < class InputIterator, class OutputIterator >
inline bool reorder_uniform( InputIterator first, InputIterator last, OutputIterator out ) {
    if( first == last ) return true;

    InputIterator tmp = first;
    double sum = 0.;
    while( tmp != last ) {
        sum -= log( *tmp );
        ++tmp;
    }

    
    double accum = 0;
    while( true ) {
        accum -= log( *first );
        ++first;

        if( first == last ) break;
        (*out) = accum / sum;
        ++out;
    }

    return true;
}

template < class InputIterator, class OutputIterator, class Order >
inline bool reorder_uniform( InputIterator first, InputIterator last, size_t N, OutputIterator out, Order & ord ) {

    double curmax = 0.;    
    while( N && first != last ) {
        curmax += (log( *first ) / ((double) N));
        (*out) = ord( exp( curmax ) );
        --N;
        ++out;
    }

    return true;
}

//template < class InputInterator, class OutputIterator >
//inline bool ascending_reorder_uniform( InputIterator first, InputIterator last, size_t N, OutputIterator out ) {
//    ascending_order ord;
//    return reorder_uniform( first, last, N, out, ord );
//}

//template < class InputIterator, class OutputIterator >
//inline bool descending_reorder_uniform( InputIterator first, InputIterator last, size_t N, OutputIterator out ) {
//    descending_order ord;
//    return reorder_uniform( first, last, N, out, ord );
//}

#endif  // REORDER_UNIFORM_HPP_
