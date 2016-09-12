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
#include "../qtl/qtl_config.hpp"

#include <iostream>

#include "qtlsim_parameters.hpp"

#include <boost/asio/io_service.hpp>
#include <boost/bind.hpp>
#include <boost/thread/thread.hpp>

#include <boost/random/mersenne_twister.hpp>

#include "qtlsim_engine_mt.hpp"

typedef boost::random::mt19937          random_engine_type;
typedef EngineMT< random_engine_type >  engine_type;

boost::mutex mio;

//typedef clotho::genetics::crossover_task< random_engine_type, unsigned long long, double > xover_type;
//
void my_task( int val ) {
    boost::lock_guard< boost::mutex > guard( mio );
    std::cout << "Val: " << val << std::endl;
}

void my_task2( int val ) {
    BOOST_LOG_TRIVIAL( info ) << "Val: " << val;
}

void run_tasks( int start, int end, int tc = 6 ) {
    boost::asio::io_service service;
    boost::thread_group threads;

    {
        std::unique_ptr< boost::asio::io_service::work > work(new boost::asio::io_service::work( service));

        // two worker threads
        for( int i = 0; i < tc; ++i )
            threads.create_thread( boost::bind( &boost::asio::io_service::run, &service) );

        for ( int i = start; i < end; ++i ) {
            service.post( boost::bind( my_task2, i ) );
        }

//        boost::this_thread::sleep_for( boost::chrono::seconds(2));
//        service.stop();
    }

    threads.join_all();
}

void writeLog( boost::property_tree::ptree & l, std::string path, std::string _type ) {
    std::shared_ptr< log_writer > tmp = makeLogWriter( path, _type );
    tmp->write( l );
}

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

    std::string path = ((out_path.empty()) ? "sample.log" : out_path + "/sample.log");
    init_logger( path, logging::trivial::info );

    boost::property_tree::ptree conf_block = config.get_child( CONFIG_BLOCK_K, config );

    generation_parameter    gen_param( conf_block );
    qtl_logging_parameter   log_param( conf_block );
    seed_parameter< >       seed_param( conf_block );

    random_engine_type  rand_engine( seed_param.m_seed );
    engine_type eng( &rand, conf_block );

    config.put_child( CONFIG_BLOCK_K, conf_block );
    writeLog( config, out_path, ".config" );

    unsigned int log_period = (( gen_param.m_size < log_param.m_period ) ? gen_param.m_size : log_param.m_period);
    simulate( eng, gen_param.m_size, log_period, out_path );

    return 0;
}

void simulate( engine_type & sim_engine, unsigned int T_gen, unsigned int log_period, std::string & out_path ) {

    std::shared_ptr< log_writer > stat_logger = makeLogWriter( out_path, ".status" );
    boost::property_tree::ptree sim_times, stat_times;

    timer_type rep_time;
    while( T_gen-- ) {
        timer_type sim_time;
        sim_engine.simulate( );
        sim_time.stop();

        clotho::utility::add_value_array( sim_times, sim_time );

        if( !( --log_period ) ) {
            boost::property_tree::ptree stat_log;
            timer_type stat_time;
            population_analyzer an( sim_engine.getChildPopulation() );
            an( stat_log );

            size_t samp_idx = 0;
            BOOST_FOREACH( auto& v, log_param.m_sampling ) {
                if( v.sample_size <= 0 ) {
                    ++samp_idx;
                    continue;
                }

                std::ostringstream oss;
                oss << "sample." << samp_idx++;

                random_sample_type samp( &rand_engine, sim_engine.getChildPopulation(), v.sample_size );

                population_analyzer _an( &samp );

                boost::property_tree::ptree samp_log;
                _an( samp_log );

                stat_log.put_child( oss.str(), samp_log );
            }

            boost::property_tree::ptree e_log;
            engine_logger_type e_logger;
            e_logger( e_log, sim_engine );

            stat_log.put_child( "population_data", e_log );

            stat_time.stop();
            clotho::utility::add_value_array( stat_times, stat_time );

            stat_logger->write( stat_log );
            log_period = ((log_param.m_period < T_gen) ? log_param.m_period : T_gen );
        }
    }

    rep_time.stop();

    boost::property_tree::ptree perform_log;
    perform_log.add( "performance.runtime", rep_time.elapsed().count() );
    perform_log.put_child( "performance.simulate", sim_times );
    perform_log.put_child( "performance.stats", stat_times );

    writeLog( perform_log, out_path, ".performance" );
}
