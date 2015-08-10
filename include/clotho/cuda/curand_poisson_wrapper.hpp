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

template < class T, class R = double >
struct fill_poisson;

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
};

}   // namespace cuda
}   // namespace clotho

#endif  // CURAND_POISSON_WRAPPER_HPP_
