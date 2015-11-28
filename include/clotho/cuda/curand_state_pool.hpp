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
#ifndef CURAND_STATE_POOL_HPP_
#define CURAND_STATE_POOL_HPP_

#include <boost/property_tree/ptree.hpp>
#include "clotho/cuda/curand_helper.hpp"

#include "clotho/random/seed_parameter.hpp"

namespace clotho {
namespace cuda {

class curand_state_pool {
public:

    typedef curandState_t       state_type;
    typedef unsigned long long  seed_type;

    typedef clotho::cuda::curand_helper< state_type > helper_type;

    static curand_state_pool * getInstance();

//    static curand_state_pool * getInstance() {
//        curand_state_pool inst;
//        return &inst;
//    }

    void initialize( boost::property_tree::ptree & config );

    state_type * get_device_states();

    unsigned int get_max_blocks() const;
    unsigned int get_max_threads() const;

    unsigned int get_total_states() const;

    virtual ~curand_state_pool();
protected:
    curand_state_pool();

    state_type  * dStates;
//    seed_type   m_seed;
    seed_parameter< seed_type > m_seed_param;

    unsigned int max_threads, max_blocks;
};

}   // namespace cuda {
}   // namespace clotho {

#endif  // CURAND_STATE_POOL_HPP_
