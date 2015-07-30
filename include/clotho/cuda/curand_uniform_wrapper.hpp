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
#ifndef CURAND_UNIFORM_WRAPPER_HPP_
#define CURAND_UNIFORM_WRAPPER_HPP_

#include <cuda.h>
#include <curand.h>
#include <thrust/device_vector.h>

#include "curand_gen_base.hpp"

namespace clotho {
namespace cuda {

template < class T >
struct fill_uniform;

template < >
struct fill_uniform < float > : public curand_gen_base {

    fill_uniform( curandGenerator_t g ) : curand_gen_base( g ) {}

    void operator()( thrust::device_vector< float > & buf ) {
        operator()( buf, buf.size() );
    }

    void operator()( thrust::device_vector< float > & buf, size_t N ) {
        float * raw_ptr = thrust::raw_pointer_cast( buf.data() );
        operator()( raw_ptr, N );
    }

    void operator()( float * buf, size_t N ) {
        if( curandGenerateUniform( gen, buf, N ) != CURAND_STATUS_SUCCESS ) {

        }
    }

    virtual ~fill_uniform() {}
};

template <>
struct fill_uniform< double > : public curand_gen_base {

    fill_uniform( curandGenerator_t g ) : curand_gen_base( g ) {}

    void operator()( thrust::device_vector< double > & buf ) {
        operator()(buf, buf.size() );
    }

    void operator()( thrust::device_vector< double > & buf, size_t N ) {
        double * raw_ptr = thrust::raw_pointer_cast( buf.data() );
        operator()( raw_ptr, N );
    }

    void operator()( double * buf, size_t N ) {
        if( curandGenerateUniformDouble( gen, buf, N ) != CURAND_STATUS_SUCCESS ) {

        }
    }

    virtual ~fill_uniform() {}
};

template <>
struct fill_uniform< unsigned int > : public curand_gen_base {
    fill_uniform( curandGenerator_t g ) : curand_gen_base( g ) {}

    void operator()( thrust::device_vector< unsigned int > & buf ) {
        operator()( buf, buf.size() );
    }

    void operator()( thrust::device_vector< unsigned int > & buf, size_t N ) {
        unsigned int * raw = thrust::raw_pointer_cast( buf.data() );
        operator()( raw, N );
    }

    void operator()( unsigned int * buf, size_t N ) {
        if( curandGenerate( gen, buf, N ) != CURAND_STATUS_SUCCESS ) {

        }
    }
};

}   // namespace cuda
}   // namespace clotho

template < class T >
void uniform_fill( curandGenerator_t g, thrust::device_vector< T > & buf, size_t N ) {
    if( buf.size() < N ) {
        buf.resize( N );
    }

    clotho::cuda::fill_uniform< T > f(g);
    f( buf, N );
}

#endif  // CURAND_UNIFORM_WRAPPER_HPP_
