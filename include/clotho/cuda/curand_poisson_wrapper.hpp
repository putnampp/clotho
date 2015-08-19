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
#ifndef CURAND_POISSON_WRAPPER_HPP_
#define CURAND_POISSON_WRAPPER_HPP_

#include <cuda.h>
#include <curand.h>
#include <thrust/device_vector.h>

#include "curand_gen_base.hpp"

namespace clotho {
namespace cuda {

template < class IntType = unsigned int, class RealType = double >
class fill_poisson : public curand_gen_base {
public:
    typedef RealType    real_type;
    typedef IntType     int_type;

    fill_poisson( curandGenerator_t g, RealType m ) : 
        curand_gen_base(g)
        , rate( m )
    {}

    inline void operator()( thrust::device_vector< IntType > & buf, size_t N ) {
        if( buf.size() < N ) { buf.resize( N ); }
        this->operator()( buf.data().get(), N );
    }

    void operator()( IntType * raw_ptr, size_t N ) {
        if( curandGeneratePoisson( this->gen, raw_ptr, N, rate ) != CURAND_STATUS_SUCCESS ) {
        }
    }

    real_type get_rate() {
        return rate;
    }

    virtual ~fill_poisson() {}

protected:
    RealType rate;
};

/*
template < >
struct fill_poisson< unsigned int, double > : public curand_gen_base {
    double mu;
    unsigned int count;

    fill_poisson( curandGenerator_t g, double m ) : curand_gen_base(g), mu( m ), count(0) {}

    void operator()( thrust::device_vector< unsigned int > & buf, size_t N ) {
//        unsigned int * raw_ptr = thrust::raw_pointer_cast( buf.data() );
        if( curandGeneratePoisson( this->gen, buf.data().get(), N, mu ) != CURAND_STATUS_SUCCESS ) {
        }
    }

    virtual ~fill_poisson() {}
};*/

}   // namespace cuda
}   // namespace clotho

#endif  // CURAND_POISSON_WRAPPER_HPP_
