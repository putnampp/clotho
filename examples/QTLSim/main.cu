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
#include "qtlsim_config.hpp"
#include <iostream>

#include "clotho/cuda/data_spaces/data_space.hpp"

#include "clotho/cuda/data_spaces/phenotype_space/device_phenotype_space.hpp"

#include "simulation_log.hpp"

#include <boost/random/mersenne_twister.hpp>

#ifdef USE_CUDA_HOST_RANDOM
#include "engine_cuda_host_random.hpp"
typedef HostPopulationSpace< real_type, int_type > population_space_type;
typedef Engine< boost::random::mt19937, real_type, int_type, int_type, cuda_host_random > engine_type;
#else
#ifdef USE_UNIT_ORDERING
typedef unit_ordered_tag< int_type > order_tag_type;
#else
typedef unordered_tag   order_tag_type;
#endif
#include "qtl_cuda_simulate_engine.hpp"
typedef PopulationSpace< real_type, int_type, order_tag_type > population_space_type;
typedef qtl_cuda_simulate_engine< population_space_type >   engine_type;

#include "clotho/cuda/data_spaces/status_buffer/status_buffer.hpp"
typedef StatusBuffer< population_space_type >       status_buffer_type;
#endif  //

static const std::string HEAP_K = "device.heap.malloc.size";

typedef std::shared_ptr< ipopulation_growth_generator >                     population_growth_generator_type;
typedef std::shared_ptr< ipopulation_growth >                               population_growth_type;

int main( int argc, char ** argv ) {

    po::variables_map vm;
    int ret = config_manager_type::getInstance()->parse_commandline( argc, argv, vm );
    if( ret ) return ret;

    boost::property_tree::ptree infile;
    getSimConfiguration( vm, infile );

    bool print_config_only = infile.empty();

    std::string prefix = "";
    if( vm.count( log_prefix_option::PREFIX_K ) ) {
        prefix = vm[ log_prefix_option::PREFIX_K ].as< std::string >();
    }

    boost::property_tree::ptree config = infile.get_child( CONFIG_BLOCK_K, infile );

    generation_parameter gen_param( config );
    unsigned int nGens = gen_param.m_size;

    std::cout << "Generations: " << nGens << std::endl;

    size_t hlimit = 0;
    cudaDeviceGetLimit( &hlimit, cudaLimitMallocHeapSize );
    hlimit = config.get< size_t >(HEAP_K, hlimit );

    assert( cudaDeviceSetLimit( cudaLimitMallocHeapSize, hlimit ) == cudaSuccess);

    config.put( HEAP_K, hlimit );

    population_growth_type pop_grow;
    if( print_config_only ) {
        population_growth_toolkit::getInstance()->tool_configurations( config );

        nGens = 0;
    } else {
        pop_grow = population_growth_toolkit::getInstance()->get_tool( config )->generate();
    }

    timer_type r;
    timer_type i_time;
    engine_type sim_engine( config );
    simulation_log log( config, prefix );
    i_time.stop();

    log.add_record( "configuration", config );

    std::string p = log.make_path( "config" );
    log.write( p );

    if( print_config_only || nGens == 0) {
        return 0;
    }

    unsigned int p_size = 0;
    unsigned int gen = 0;

    boost::property_tree::ptree _sim, _an, _an_state, _an_samps;
    boost::property_tree::ptree free_count, a_count, var_count;

#ifndef USE_CUDA_HOST_RANDOM
    status_buffer_type stats( nGens );
#endif  // USE_CUDA_HOST_RANDOM

    while( gen < nGens ) {
        assert( gen < 10);
        p_size = (*pop_grow)( p_size, gen++ );

        timer_type t;
        sim_engine.simulate(p_size);
        t.stop();

#ifdef USE_CUDA_HOST_RANDOM
        sim_engine.trackStats();
#else
        stats.track( sim_engine.get_offspring_population(), (gen - 1) );
#endif  // USE_CUDA_HOST_RANDOM

        clotho::utility::add_value_array( _sim, t );

        if( log.hasPeriodElapsed() ) {

            // population level analysis
            t.start();
            sim_engine.analyze_population();
            t.stop();

            clotho::utility::add_value_array( _an, t );

            t.start();
            log.record_state( &sim_engine );
            t.stop();

            clotho::utility::add_value_array( _an_state, t );

            typedef std::vector< sample_log_params > samples_type;
            typedef typename samples_type::iterator  sample_iter;

            t.start();
            boost::property_tree::ptree samps;
            unsigned int i = 0;
            for( sample_iter it = log.m_sampling.begin(); it != log.m_sampling.end(); it++ ) {
                if( it->pairwise ) {
                    boost::property_tree::ptree res;
                    sim_engine.analyze_sample( it->sample_size, res );

                    std::ostringstream oss;
                    oss << i++;
                    samps.add_child( oss.str(), res );
                }
            }
            t.stop();

            clotho::utility::add_value_array( _an_samps, t);

            if(!samps.empty() ) {
                log.add_record( "samples", samps );
            }

            //unsigned int fc = log.getLog().get< unsigned int>( "population.current.free_space.total", 0);
            //unsigned int ac = log.getLog().get< unsigned int>( "population.current.alleles.capacity", 0);

            //clotho::utility::add_value_array(free_count, fc);
            //clotho::utility::add_value_array(a_count, ac );
            //clotho::utility::add_value_array(var_count, (ac - fc));
            
//            std::cerr << gen << ": Writing log file" << std::endl;
            log.write();
        }
    }
    r.stop();

    boost::property_tree::ptree _run, _mem;
    _run.put( "total", r );
    _run.put( "init", i_time );
    _run.put_child( "simulate", _sim );
    _run.put_child( "analysis.global.elapsed", _an );
    _run.put_child( "analysis.global.record_elapsed", _an_state );
    _run.put_child( "analysis.samples.elapsed", _an_samps );

    log.add_record( "performance", _run );

//    _mem.put_child( "free_count", free_count);
//    _mem.put_child( "allele_count", a_count);
//    _mem.put_child( "variable_count", var_count);

#ifdef USE_CUDA_HOST_RANDOM
    sim_engine.get_tracked_states( _mem );
#else
    stats.get_state( _mem );
#endif // USE_CUDA_HOST_RANDOM
    log.add_record( "memory", _mem );
    p = log.make_path( "performance" );
    log.write( p );

    return 0;
}
