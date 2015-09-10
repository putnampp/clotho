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
#ifndef QTL_CUDA_SIMULATE_ENGINE_HPP_
#define QTL_CUDA_SIMULATE_ENGINE_HPP_

#include <algorithm>
#include <boost/property_tree/ptree.hpp>

#include "simulation_log.hpp"

//#include "clotho/cuda/mutation/population_mutation_generator.hpp"
#include "clotho/cuda/data_spaces/event_space/device_event_space.hpp"
#include "clotho/cuda/mutation/mutation_event_generator.hpp"
#include "clotho/cuda/crossover/crossover_generator.hpp"
#include "clotho/cuda/selection/selection_event_generator.hpp"

#include "clotho/utility/timer.hpp"

template < class PopulationSpaceType > 
class qtl_cuda_simulate_engine {
public:
    typedef PopulationSpaceType                                 population_space_type;

    typedef typename population_space_type::allele_space_type   allele_space_type;
    typedef typename population_space_type::sequence_space_type sequence_space_type;

//    typedef population_mutation_generator< population_space_type >  mutation_generator_type;
//    typedef MutationGenType                                     mutation_generator_type;
//
//
    typedef typename allele_space_type::real_type                       real_type;
    typedef typename allele_space_type::int_type                        int_type;
    typedef typename allele_space_type::order_tag_type                  order_tag_type;

    typedef MutationEventGenerator< real_type, int_type, order_tag_type >   mutation_generator_type;
    typedef typename mutation_generator_type::space_type            mutation_event_space_type;

    typedef CrossoverGenerator< population_space_type >             crossover_generator_type;

    typedef SelectionEventGenerator                                 selection_generator_type;

    typedef clotho::utility::timer                                  timer_type;

    qtl_cuda_simulate_engine( boost::property_tree::ptree & config ) :
        current_pop( &hPop0 )
        , prev_pop( &hPop1 )
        , m_log( config )
        , mut_gen( config )
        , xover_gen( config )
        , sel_gen(config)
    {
        parse_config( config );
        initialize();
    }

    void simulate( unsigned int N ) {

        timer_type t;
        unsigned int cur_seq_count = 2 * N;

        mut_gen.generate( dMutations, cur_seq_count );

//        cudaDeviceSynchronize();

//        std::cerr << "Events: " << dMutations << std::endl;

        current_pop->resize( prev_pop, dMutations, cur_seq_count );

        xover_gen( current_pop );

        // recombine parents 
        sel_gen.generate( prev_pop, current_pop );

        mut_gen.scatter( current_pop, dMutations, cur_seq_count );

        cudaDeviceSynchronize();
        t.stop();
//        std::cerr << "Current Population: " << *current_pop << std::endl;
        std::cerr << "Lapse: " << t << std::endl;
        swap();
    }

    void swap() {
        std::swap( prev_pop, current_pop );
    }

    virtual ~qtl_cuda_simulate_engine() {
        delete_space( dMutations );
    }

protected:

    void parse_config( boost::property_tree::ptree & config ) {
    }

    void initialize( ) {
        create_space( dMutations );
    }

    mutation_event_space_type     * dMutations;
    population_space_type   hPop0, hPop1;
    population_space_type   * prev_pop, * current_pop;

    simulation_log              m_log;
    mutation_generator_type     mut_gen;
    crossover_generator_type    xover_gen;
    selection_generator_type    sel_gen;
};

#endif  // QTL_CUDA_SIMULATE_ENGINE_HPP_
