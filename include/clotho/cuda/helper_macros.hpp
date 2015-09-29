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
#ifndef CLOTHO_CUDA_HELPER_MACROS_HPP_
#define CLOTHO_CUDA_HELPER_MACROS_HPP_

#ifdef CUDA_CLOTHO_SAFE_KERNELS

#define CHECK_LAST_KERNEL_EXEC                                                  \
{                                                                               \
    cudaDeviceSynchronize();                                                    \
    cudaError_t err = cudaGetLastError();                                       \
    if( err != cudaSuccess ) {                                                  \
        std::cerr << "CUDA error: " << cudaGetErrorString( err ) << std::endl;  \
        assert(false);                                                          \
    }                                                                           \
}

#else

#define CHECK_LAST_KERNEL_EXEC

#endif  // CUDA_CLOTHO_SAFE_KERNELS

#endif  // CLOTHO_CUDA_HELPER_MACROS_HPP_
