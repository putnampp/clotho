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
#ifndef SAMPLE_POPULATION_HPP_
#define SAMPLE_POPULATION_HPP_

#include "clotho/cuda/data_spaces/population_space/population_space_def.hpp"
#include "clotho/cuda/data_spaces/basic_data_space.hpp"

#include "clotho/cuda/curand_state_pool.hpp"
#include "clotho/cuda/sampling/random_sample.hpp"

template < class PopulationSpaceType >
struct SamplePopulation;

template < class RealType, class IntType, class OrderTag >
struct SamplePopulation< PopulationSpace< RealType, IntType, OrderTag > > {
    typedef PopulationSpace< RealType, IntType, OrderTag >  population_type;
    typedef clotho::cuda::curand_state_pool                 state_pool_type;

    population_type                     * m_pop;
    basic_data_space< unsigned int >    * m_dSubpop;
    SamplePopulation( population_type * p, unsigned int N = 0 ) :
        m_pop( NULL )
    {
        create_space( m_dSubpop );

        set_population( p, N );
    }

    void set_population( population_type * p, unsigned int N ) {
        m_pop = p;

        N = ((m_pop) ? N : 0 );

        random_sample<<< 1, 1024 >>>( state_pool_type::getInstance()->get_device_states(), m_pop->sequences.get_device_space(), N, m_dSubpop );
        CHECK_LAST_KERNEL_EXEC
    }

    virtual ~SamplePopulation() {
        delete_space( m_dSubpop );
    }
};

#endif  // SAMPLE_POPULATION_HPP_
