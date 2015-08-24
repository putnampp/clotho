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
#ifndef POISSON_DISTRIBUTION_HPP_
#define POISSON_DISTRIBUTION_HPP_

#include <cuda.h>

template < class RealType, unsigned int K = 32 >
struct poisson_cdf {
    static const unsigned int MAX_K = K;
    RealType _cdf[ MAX_K ];
    unsigned int max_k;
};

template < class RealType >
__global__ void  make_poisson_cdf_maxk32( poisson_cdf< RealType, 32 > * c, RealType lambda ) {
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int lane_id = (tid & 31);

    RealType p = ((lane_id != 0) ? (lambda / (RealType) lane_id) : 1.0);

    for( unsigned int i = 1; i < 32; i <<= 1 ) {
        RealType _c = __shfl_up(p, i );
        p *=  ((lane_id >= i) ? _c : 1.0);
    }

    for( unsigned int i = 1; i < 32; i <<= 1 ) {
        RealType _c = __shfl_up( p, i );
        p += ((lane_id >= i ) ? _c : 0.0 );
    }

    p *= exp( 0.0 - lambda ); // $p_{lane_id} = \sum_{k=0}^{lane_id} \frac{\lambda^{k}}{lane_id!} e^{-\lambda} $

    c->_cdf[ tid ] = p;

    unsigned int koff = tid;
    for( unsigned int i = 1; i < 32; i <<= 1 ) {
        RealType _c = __shfl_up( p, i );
        unsigned int _k = __shfl_up( koff, i);

        koff = (((_c >= 1.0) && (_k < koff)) ? _k : koff);
    }

    if( tid == 31 ) {
        c->max_k = koff;
    }
}

template < class RealType >
__host__ __device__ unsigned int _find_poisson_maxk32( volatile RealType * cdf, RealType x, unsigned int N = 32 ) {
    unsigned int k = 0;
    while( N-- ) {
        RealType _c = cdf[k];
        //k += (_c < x);
        k += (( _c < x) ? 1 : 0);
    }
    return k;
}
#endif  // POISSON_DISTRIBUTION_HPP_