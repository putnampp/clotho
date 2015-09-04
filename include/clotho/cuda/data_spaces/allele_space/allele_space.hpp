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

#include "clotho/cuda/data_spaces/allele_space/device_allele_space.hpp"

template < class RealType, class IntType, class OrderTag >
class AlleleSpace {
public:
    typedef device_allele_space< RealType, IntType, OrderTag >  device_space_type;

    typedef typename device_space_type::real_type   real_type;
    typedef typename device_space_type::int_type    int_type;
    typedef typename device_space_type::order_tag_type  order_tag_type;

    typedef AlleleSpace< RealType, IntType, OrderTag > self_type;

    AlleleSpace( );

    unsigned int total_free_space();

    void resize( unsigned int N );

    device_space_type * get_device_space();

//    void merge_alleles( self_type & alls, self_type * muts);

    void expand_relative_to( self_type & alls, device_event_space< IntType, OrderTag > * space );

    virtual ~AlleleSpace();

    template < class R, class I, class O >
    friend std::ostream & operator<<( std::ostream & out, const AlleleSpace< R, I, O> & s );

protected:
    void initialize();

    device_space_type  * dAlleles;    // exists on device
};


#define _HEADER template < class RealType, class IntType, class OrderTag >
#define _CLASS  AlleleSpace< RealType, IntType, OrderTag >

_HEADER
_CLASS::AlleleSpace() : dAlleles( NULL ) {
    initialize();
}

_HEADER
void _CLASS::initialize() {
    create_space( dAlleles, 0 );
}

_HEADER
unsigned int _CLASS::total_free_space() {
    return 0;
}

_HEADER
void _CLASS::resize( unsigned int N ) {
    resize_space( dAlleles, N );
}

_HEADER
typename _CLASS::device_space_type * _CLASS::get_device_space() {
    return dAlleles;
}

/*
_HEADER
void _CLASS::merge_alleles( self_type * A, self_type * B ) {
    merge_allele_space( A->dAlleles, B->dAlleles, dAlleles );
}*/

_HEADER
void _CLASS::expand_relative_to( self_type & alls, device_event_space< IntType, OrderTag > * muts ) {
    merge_space( alls.dAlleles, muts, dAlleles );
}

_HEADER
_CLASS::~AlleleSpace( ) {
    delete_space( dAlleles );
}

_HEADER
std::ostream & operator<<( std::ostream & out, const _CLASS & rhs ) {
    typedef typename _CLASS::device_space_type space_type;
    space_type hCopy;

    assert( cudaMemcpy( &hCopy, rhs.dAlleles, sizeof( space_type ), cudaMemcpyDeviceToHost ) == cudaSuccess );

    out << hCopy;

    return out;
}

#undef _CLASS
#undef _HEADER

#endif  // ALLELE_SPACE_HPP_
