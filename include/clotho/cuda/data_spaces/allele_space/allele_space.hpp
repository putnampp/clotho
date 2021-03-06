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

#include <boost/property_tree/ptree.hpp>
#include "clotho/utility/clotho_strings.hpp"

#include "clotho/cuda/data_spaces/allele_space/device_allele_space.hpp"
#include "clotho/cuda/data_spaces/free_space/device_free_space.hpp"
#include "clotho/cuda/data_spaces/event_space/device_event_space.hpp"
#include "clotho/cuda/device_state_object.hpp"
#include "clotho/cuda/helper_macros.hpp"

template < class RealType >
class AlleleSpace : public clotho::utility::iStateObject {
public:
//    typedef device_allele_space< RealType >  device_space_type;

    typedef device_weighted_allele_space< RealType >    device_space_type;

    typedef typename device_space_type::real_type   real_type;

    typedef AlleleSpace< RealType > self_type;

    AlleleSpace( boost::property_tree::ptree & config );

    void setNeutralityP( real_type p );

    unsigned int total_free_space();

    void resize( unsigned int N );

    device_space_type * get_device_space();

    template < class IntType, class OrderTag >
    void expand_relative_to( self_type & alls, device_free_space< IntType, OrderTag > * fspace, device_free_space< IntType, OrderTag > * ofspace, device_event_space< IntType, OrderTag > * space );

    void get_state( boost::property_tree::ptree & state );

    template < class IntType, class OrderTag >
    void move_fixed( self_type & alls, device_free_space< IntType, OrderTag > * space );

    virtual ~AlleleSpace();

    template < class R >
    friend std::ostream & operator<<( std::ostream & out, const AlleleSpace< R > & s );

protected:
    void initialize();

    device_space_type  * dAlleles;    // exists on device
};


//#define _HEADER template < class RealType, class IntType, class OrderTag >
//#define _CLASS  AlleleSpace< RealType, IntType, OrderTag >
//
#define _HEADER template < class RealType >
#define _CLASS  AlleleSpace< RealType >

_HEADER
_CLASS::AlleleSpace( boost::property_tree::ptree & config ) : dAlleles( NULL ) {

    boost::property_tree::ptree /*all,*/ neu;

//    all = config.get_child( ALLELE_BLOCK_K, all );
    neu = config.get_child(NEUTRAL_BLOCK_K, neu );

    real_type p = neu.get< real_type >( P_K, 0.0 );

    assert( 0.0 <= p && p <= 1.0 );

    neu.put( P_K, p );
    config.put_child( NEUTRAL_BLOCK_K, neu );
//    config.put_child( ALLELE_BLOCK_K, neu );

    initialize();

    setNeutralityP( p );
}

_HEADER
void _CLASS::initialize() {
    create_space( dAlleles );
}

_HEADER
void _CLASS::setNeutralityP( real_type p ) {
    typedef typename _CLASS::device_space_type space_type;
    space_type hCopy;

    assert( cudaMemcpy( &hCopy, dAlleles, sizeof( space_type ), cudaMemcpyDeviceToHost ) == cudaSuccess );

    hCopy.neutral_p = p;
    assert( cudaMemcpy( dAlleles, &hCopy, sizeof( space_type ), cudaMemcpyHostToDevice ) == cudaSuccess );
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

_HEADER template < class IntType, class OrderTag >
void _CLASS::expand_relative_to( self_type & alls, device_free_space< IntType, OrderTag > * fspace, device_free_space< IntType, OrderTag > * ofspace, device_event_space< IntType, OrderTag > * muts ) {
    merge_space( alls.dAlleles, fspace, muts, ofspace, dAlleles );
}

_HEADER template < class IntType, class OrderTag >
void _CLASS::move_fixed( self_type & alls, device_free_space< IntType, OrderTag > * free_space ) {

    move_fixed_allele_kernel<<< 1, 32 >>>( alls.dAlleles, dAlleles, free_space );
    CHECK_LAST_KERNEL_EXEC
}

_HEADER
void _CLASS::get_state( boost::property_tree::ptree & state ) {
    get_device_object_state( state, dAlleles );
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
