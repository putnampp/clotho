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
#ifndef DATA_SPACE_HELPER_HPP_
#define DATA_SPACE_HELPER_HPP_

template < class DataType >
__global__ void copy_heap( DataType * d_loc, DataType * d_heap, unsigned int N ) {
    memcpy( d_loc, d_heap, N * sizeof( DataType ) );
}

template < class DataType >
void copy_heap_data( DataType * hMem, DataType * dMem, unsigned int N ) {
    DataType * tmp;

    size_t sBytes = N * sizeof( DataType );
    assert( cudaMalloc( (void **) &tmp, sBytes ) == cudaSuccess );

    copy_heap<<< 1, 1>>>( tmp, dMem, N );
    cudaDeviceSynchronize();

    assert( cudaMemcpy( hMem, tmp, sBytes, cudaMemcpyDeviceToHost ) == cudaSuccess );

    cudaFree( tmp );
}

#endif  // DATA_SPACE_HELPER_HPP_
