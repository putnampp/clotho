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
#ifndef CURAND_GEN_BASE_HPP_
#define CURAND_GEN_BASE_HPP_

#include <curand.h>

namespace clotho {
namespace cuda {

struct curand_gen_base {
    curandGenerator_t gen;

    curand_gen_base( curandGenerator_t g ) : gen(g) {}

    virtual ~curand_gen_base() {}
};

}   // namespace cuda
}   // namespace clotho

#endif  // CURAND_GEN_BASE_HPP_
