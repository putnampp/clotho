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

#include <iostream>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "clotho/configuration_manager/config_manager.hpp"

#include "clotho/cuda/data_spaces/data_space.hpp"

#include "clotho/cuda/data_spaces/phenotype_space/device_phenotype_space.hpp"

#include "options/configuration_option.hpp"
#include "options/log_prefix_option.hpp"

#include "qtl_cuda_simulate_engine.hpp"
#include "clotho/genetics/population_growth_toolkit.hpp"

#include "simulation_log.hpp"

#include "clotho/utility/timer.hpp"
#include "clotho/utility/log_helper.hpp"

typedef clotho::configuration_manager::config_manager       config_manager_type;

#ifdef USE_DOUBLE_REAL
typedef double          real_type;
#else
typedef float           real_type;
#endif
typedef unsigned int    int_type;

#ifdef USE_UNIT_ORDERING
typedef unit_ordered_tag< int_type > order_tag_type;
#else
typedef unordered_tag   order_tag_type;
#endif

typedef PopulationSpace< real_type, int_type, order_tag_type > population_space_type;

typedef qtl_cuda_simulate_engine< population_space_type >   engine_type;

typedef clotho::utility::timer                              timer_type;

static const std::string GEN_K = "generations";
static const std::string HEAP_K = "device.heap.malloc.size";

typedef std::shared_ptr< ipopulation_growth_generator >                     population_growth_generator_type;
typedef std::shared_ptr< ipopulation_growth >                               population_growth_type;

int main( int argc, char ** argv ) {

    po::variables_map vm;
    int ret = config_manager_type::getInstance()->parse_commandline( argc, argv, vm );
    if( ret ) return ret;

    boost::property_tree::ptree infile;
    bool print_config_only = true;
    if( vm.count( configuration_option::CONFIG_K ) ) {
        std::string p = vm[ configuration_option::CONFIG_K ].as< configuration_option::path_type >();

        boost::property_tree::read_json( p, infile );
        print_config_only = false;
    }

    std::string prefix;
    if( vm.count( log_prefix_option::PREFIX_K ) ) {
        prefix = vm[ log_prefix_option::PREFIX_K ].as< std::string >();
    }

    boost::property_tree::ptree config = infile.get_child( "configuration", infile );

    unsigned int nGens = 100;

    if( config.get_child_optional( GEN_K ) == boost::none ) {
        config.put(GEN_K, nGens);
    } else {
        nGens = config.get< unsigned int >( GEN_K, nGens);
    }

    if( config.get_child_optional( HEAP_K ) == boost::none ) {
        size_t hlimit = 0;
        cudaDeviceGetLimit( &hlimit, cudaLimitMallocHeapSize );

        config.put( HEAP_K, hlimit );
    } else {
        size_t hlimit = config.get< size_t >( HEAP_K );

        cudaDeviceSetLimit( cudaLimitMallocHeapSize, hlimit );
    }

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
    simulation_log log( config );
    i_time.stop();

    if( !prefix.empty() ) {
        log.set_path_prefix( prefix );
    }

    log.add_record( "configuration", config );

    std::string p = log.make_path( "config" );
    log.write( p );

    if( print_config_only ) {
        return 0;
    }

    unsigned int p_size = 0;
    unsigned int gen = 0;

    boost::property_tree::ptree _sim, _an;
    while( gen < nGens ) {
        p_size = (*pop_grow)( p_size, gen++ );

        timer_type t;
        sim_engine.simulate(p_size);
        t.stop();

        clotho::utility::add_value_array( _sim, t );

        t.start();
        sim_engine.analyze_population();
        t.stop();

        clotho::utility::add_value_array( _an, t );

        if( log( &sim_engine ) ) {
            std::cerr << gen << ": Writing log file" << std::endl;
            log.write();
        }
    }
    r.stop();

    boost::property_tree::ptree _run;
    _run.put( "total", r );
    _run.put( "init", i_time );
    _run.put_child( "simulate", _sim );
    _run.put_child( "analyze", _an );

    log.add_record( "runtime", _run );

    p = log.make_path( "profile" );
    log.write( p );

    return 0;
}
