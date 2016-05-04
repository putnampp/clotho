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

#include "log_writer.hpp"

#include <iostream>
#include <algorithm>
#include <map>

#include <boost/random/mersenne_twister.hpp>

#include "qtlsim_engine.hpp"
#include "clotho/random/seed_parameter.hpp"
#include "../qtl/qtl_logging_parameter.hpp"

typedef boost::random::mt19937          random_engine_type;
typedef Engine< random_engine_type >    simulate_engine_type;


//#define XSTR( x ) #x
//#define STR( x )  XSTR( x )
//
//const std::string ENGINE_BLOCK_K = "engine";
//const std::string POWERSET_BLOCK_K = "powerset";
//const std::string SUBSET_BLOCK_K = "subset";
//const std::string CLASSIFIER_BLOCK_K = "classifier";
//void write_engine_config( boost::property_tree::ptree & elog );

int main( int argc, char ** argv ) {

    po::variables_map vm;
    int ret = config_manager_type::getInstance()->parse_commandline( argc, argv, vm );
    if( ret ) return ret;

    random_engine_type  rand;
    boost::property_tree::ptree config;

    getSimConfiguration( vm, config );

    bool print_config_only = config.empty();

    std::string out_path = "";
    if( vm.count( log_prefix_option::PREFIX_K ) ) {
        out_path = vm[ log_prefix_option::PREFIX_K ].as< std::string >();
    }

    std::shared_ptr< log_writer > config_logger = makeLogWriter( out_path, ".config" );
    std::shared_ptr< log_writer > performance_logger = makeLogWriter( out_path, ".performance" );
    std::shared_ptr< log_writer > stat_logger = makeLogWriter( out_path, ".status" );

    boost::property_tree::ptree conf_block = config.get_child( CONFIG_BLOCK_K, config );

    generation_parameter    gen_param( conf_block );
    qtl_logging_parameter   log_param( conf_block );
    seed_parameter< >       seed_param( conf_block );

    random_engine_type  rand_engine( seed_param.m_seed );

    simulate_engine_type sim_engine( &rand_engine, conf_block );

    config.put_child( CONFIG_BLOCK_K, conf_block );
    config_logger->write( config );

    if( !print_config_only ) {
        unsigned int log_period = (( gen_param.m_size < log_param.m_period ) ? gen_param.m_size : log_param.m_period);

        boost::property_tree::ptree sim_times, stat_times;

        timer_type rep_time;
        unsigned int T_gen = gen_param.m_size;
        while( --T_gen ) {
            timer_type sim_time;
            //sim_engine.simulate( );
            sim_time.stop();

            clotho::utility::add_value_array( sim_times, sim_time );

            if( !( --log_period ) ) {
                boost::property_tree::ptree stat_log;
                timer_type stat_time;

                log_period = ((log_param.m_period < T_gen) ? log_param.m_period : T_gen - 1 );
                stat_time.stop();

                clotho::utility::add_value_array( stat_times, stat_time );

                stat_logger->write( stat_log );
            }
        }

        rep_time.stop();

        boost::property_tree::ptree perform_log;
        perform_log.add( "performance.runtime", rep_time.elapsed().count() );
        perform_log.put_child( "performance.simulate", sim_times );
        perform_log.put_child( "performance.stats", stat_times );

        performance_logger->write( perform_log );
    }

    return 0;
}

//void write_engine_config( boost::property_tree::ptree & elog ) {
//    boost::property_tree::ptree sset, recomb, pset;
//    recomb.put( "tag0", STR( RECOMBINE_INSPECT_METHOD ) );
//    recomb.put( "tag1", STR( BIT_WALK_METHOD ) );
//    sset.put_child( "recombine", recomb );
//    sset.put( TYPE_K, STR( SUBSETTYPE ) );
//
//    pset.put_child( SUBSET_BLOCK_K, sset );
////    pset.put( "block_size", BLOCK_UNIT_SIZE );
////    compile_log.put(ENGINE_BLOCK_K + "." + REC_BLOCK_K + "." + TYPE_K, STR( (RECOMBTYPE) ) );
////    compile_log.put(ENGINE_BLOCK_K + "." + POWERSET_BLOCK_K +"." + SIZE_K, BLOCK_UNIT_SIZE );
//
////    compile_log.put(ENGINE_BLOCK_K + "." + POWERSET_BLOCK_K + "." + SUBSET_BLOCK_K  +"." + TYPE_K, STR( SUBSETTYPE ) );
////    compile_log.put( ENGINE_BLOCK_K + "." + POWERSET_BLOCK_K + "." + SUBSET_BLOCK_K + "." + REC_BLOCK_K + ".tag0", STR( RECOMBINE_INSPECT_METHOD ) );
////    compile_log.put( ENGINE_BLOCK_K + "." + POWERSET_BLOCK_K + "." + SUBSET_BLOCK_K + "." + REC_BLOCK_K + ".tag1", STR( BIT_WALK_METHOD ) );
//
////    compile_log.put( ENGINE_BLOCK_K + ".reproduction_method.type", STR(REPRODUCTION_METHOD_TAG));
////
////    compile_log.put( ENGINE_BLOCK_K + ".individual_selector.type", STR(IND_SELECT) );
//
//    boost::property_tree::ptree eng;
//    eng.put_child( POWERSET_BLOCK_K, pset );
//    eng.put( REC_BLOCK_K + ".type", STR( (RECOMBTYPE) ) );
//    eng.put( "reproduction_method.type", STR(REPRODUCTION_METHOD_TAG));
//    eng.put( "individual_selector.type", STR(IND_SELECT) );
//    eng.put( "description", "Simulator compiled objects; READ ONLY");
//
//    elog.put_child( ENGINE_BLOCK_K, eng );
//}
