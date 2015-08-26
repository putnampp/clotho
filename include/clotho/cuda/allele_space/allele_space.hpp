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
#ifndef ALLELE_SPACE_HPP_
#define ALLELE_SPACE_HPP_

#include "clotho/cuda/allele_space/device_allele_space.hpp"

template < class RealType, class IntType, class OrderTag >
class AlleleSpace {
public:
    typedef device_allele_space< RealType, IntType, OrderTag > device_space_type;

    typedef AlleleSpace< RealType, IntType, OrderTag > self_type;

    AlleleSpace( );

    unsigned int total_free_space();

    void resize( unsigned int N );

    virtual ~AlleleSpace();

    template < class R, class I, class O >
    friend std::ostream & operator<<( std::ostream & out, const AlleleSpace< R, I, O> & s );

protected:
    void initialize();

    device_space_type  * dSpace;    // exists on device
    device_space_type hSpace;       // exists on host
};


#define _HEADER template < class RealType, class IntType, class OrderTag >
#define _CLASS  AlleleSpace< RealType, IntType, OrderTag >

_HEADER
_CLASS::AlleleSpace() : dSpace( NULL ) {
    initialize();
}

_HEADER
void _CLASS::initialize() {
    allele_space_alloc( dSpace, 0 );
}

_HEADER
unsigned int _CLASS::total_free_space() {
    return 0;
}

_HEADER
void _CLASS::resize( unsigned int N ) {
    resize_space( dSpace, N );
}

_HEADER
_CLASS::~AlleleSpace( ) {
    allele_space_free( dSpace );
}

_HEADER
std::ostream & operator<<( std::ostream & out, const _CLASS & rhs ) {
    typedef typename _CLASS::device_space_type space_type;
    space_type hCopy;

    assert( cudaMemcpy( &hCopy, rhs.dSpace, sizeof( space_type ), cudaMemcpyDeviceToHost ) == cudaSuccess );

    out << hCopy;

    return out;
}

#undef _CLASS
#undef _HEADER

#endif  // ALLELE_SPACE_HPP_
