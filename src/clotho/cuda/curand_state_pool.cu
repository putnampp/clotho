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
#include "clotho/cuda/curand_state_pool.hpp"

#include "clotho/utility/timer.hpp"
#include <iostream>

namespace clotho {
namespace cuda {

curand_state_pool::curand_state_pool() :
    dStates( NULL )
    , m_seed_param(0)
    , max_threads(32)
    , max_blocks(1)
{}

curand_state_pool * curand_state_pool::getInstance() {
    static curand_state_pool inst;
    return &inst;
}

void curand_state_pool::initialize( boost::property_tree::ptree & config ) {
    if( dStates != NULL ) return;

    seed_parameter< seed_type > tmp_seed(config);

    m_seed_param = tmp_seed;

    boost::property_tree::ptree lconfig;

    if( config.get_child_optional( "state_pool" ) != boost::none ) {
        lconfig = config.get_child( "state_pool" );
    }

/*
    if( lconfig.get_child_optional( "seed" ) != boost::none ) {
        m_seed = lconfig.get< seed_type >( "seed", m_seed );
    }

    if( m_seed == 0 ) {
        m_seed = clotho::utility::clock_type::now().time_since_epoch().count();
        lconfig.put("seed", m_seed );
    }
*/

    if( lconfig.get_child_optional( "max_blocks" ) == boost::none ) {
        lconfig.put("max_blocks", max_blocks );
    } else {
        max_blocks = lconfig.get< unsigned int >( "max_blocks", max_blocks );
    }

    if( lconfig.get_child_optional( "max_threads" ) == boost::none) {
        lconfig.put("max_threads", max_threads );
    } else {
        max_threads = lconfig.get< unsigned int >("max_threads", max_threads );
    }

    config.put_child( "state_pool", lconfig );

    if( max_threads * max_blocks != 0 ) {
        std::cerr << "Seed Value: " << m_seed_param.m_seed << std::endl;
        std::cerr << "Initializing state pool: " << max_blocks << ", " << max_threads << std::endl;
        helper_type::make_states( dStates, m_seed_param.m_seed, max_blocks, max_threads );
    }
}

typename curand_state_pool::state_type * curand_state_pool::get_device_states() {
    return dStates;
}

unsigned int curand_state_pool::get_max_blocks() const {
    return max_blocks;
}

unsigned int curand_state_pool::get_max_threads() const {
    return max_threads;
}

unsigned int curand_state_pool::get_total_states() const {
    return max_threads * max_blocks;
}

curand_state_pool::~curand_state_pool() {
    helper_type::cleanup_states( dStates );
    dStates = NULL;
}

}   // namespace cuda {
}   // namespace clotho {
