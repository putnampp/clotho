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
#ifndef INDEX_CONVERTER_HPP_
#define INDEX_CONVERTER_HPP_

#include <cuda.h>

template < unsigned int WIDTH >
struct fixed_width_converter {

    __forceinline__ __device__ unsigned int major_offset( unsigned int idx ) {
        return idx / WIDTH;
    }

    __forceinline__ __device__ unsigned int minor_offset( unsigned int idx ) {
        return idx % WIDTH;
    }
};

/*
template <  >
struct fixed_width_converter< 64 > {

    __forceinline__ __device__ unsigned int major_offset( unsigned int idx ) {
        return idx >> 6;
    }

    __forceinline__ __device__ unsigned int minor_offset( unsigned int idx ) {
        return idx & 63;
    }
};*/

template <  >
struct fixed_width_converter< 32 > {

    __forceinline__ __device__ unsigned int major_offset( unsigned int idx ) {
        return (idx >> 5);
    }

    __forceinline__ __device__ unsigned int minor_offset( unsigned int idx ) {
        return (idx & 31);
    }
};

/*
template <  >
struct fixed_width_converter< 16 > {

    __forceinline__ __device__ unsigned int major_offset( unsigned int idx ) {
        return idx >> 4;
    }

    __forceinline__ __device__ unsigned int minor_offset( unsigned int idx ) {
        return idx & 15;
    }
};

template <  >
struct fixed_width_converter< 8 > {

    __forceinline__ __device__ unsigned int major_offset( unsigned int idx ) {
        return idx >> 3;
    }

    __forceinline__ __device__ unsigned int minor_offset( unsigned int idx ) {
        return idx & 7;
    }
};*/
#endif  // INDEX_CONVERTER_HPP_
