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
#include <set>
#include <sstream>

#include <boost/random/mersenne_twister.hpp>

#include "engines.hpp"

#include "qtlsim_parameters.hpp"

#include "clotho/data_spaces/generators/random_sample_generator.hpp"
#include "clotho/data_spaces/analysis/allele_frequency/allele_frequency.hpp"
#include "clotho/data_spaces/analysis/pairwise_difference/pairwise_difference.hpp"

typedef boost::random::mt19937          random_engine_type;

typedef Engine< random_engine_type, real_type, block_unit_type, unsigned int, ENGINE_MODE >    simulate_engine_type;
typedef typename simulate_engine_type::sequence_space_type      sequence_space_type;
typedef typename simulate_engine_type::allele_type              allele_space_type;

typedef clotho::genetics::random_sample_generator< random_engine_type >  random_sample_type;
typedef clotho::genetics::allele_frequency< sequence_space_type, unsigned int > allele_freq_type;
typedef clotho::genetics::pairwise_difference< sequence_space_type > pairwise_diff_type;

typedef clotho::utility::state_getter< simulate_engine_type >   engine_logger_type;

#ifdef DEBUG_MODE
#define CLOTHO_LOG_LEVEL logging::trivial::debug
#else
#define CLOTHO_LOG_LEVEL logging::trivial::info
#endif  // DEBUG_MODE

/**
 * Population Analyzer
 *
 * Helper object to group the population statistics to be periodically evaluated
 *
 */
class population_analyzer {
public:
    population_analyzer( random_sample_type * s, sequence_space_type * ss, allele_space_type * alls );

    void operator()( boost::property_tree::ptree & log );

    virtual ~population_analyzer() {}

protected:
    sequence_space_type * m_seqs;
    allele_space_type   * m_alleles;
    random_sample_type * m_samp;
};

void writeLog( boost::property_tree::ptree & l, std::string path, std::string _type ) {
    std::shared_ptr< log_writer > tmp = makeLogWriter( path, _type );
    tmp->write( l );
}


/**
 * QTLSimMT
 *
 * An example application built using the Clotho API.  
 *
 * The main application performs the following steps:
 *
 *   - Parse the input commandline arguments; use -h to retrieve a listing of the available options
 *   - Setup basic simulation configuration (ie. logging paths; simulation engine;)
 *   - 
 */
int main( int argc, char ** argv ) {

    po::variables_map vm;
    int ret = config_manager_type::getInstance()->parse_commandline( argc, argv, vm );
    if( ret ) return ret;

    boost::property_tree::ptree config;

    getSimConfiguration( vm, config );

    bool print_config_only = config.empty();

    std::string out_path = "";
    if( vm.count( log_prefix_option::PREFIX_K ) ) {
        out_path = vm[ log_prefix_option::PREFIX_K ].as< std::string >();
    }

    std::string log_path = out_path + ".log";
    init_logger( log_path, CLOTHO_LOG_LEVEL );

    std::shared_ptr< log_writer > stat_logger = makeLogWriter( out_path, ".status" );

    boost::property_tree::ptree conf_block = config.get_child( CONFIG_BLOCK_K, config );

    generation_parameter    gen_param( conf_block );
    qtl_logging_parameter   log_param( conf_block );
    seed_parameter< >       seed_param( conf_block );

    random_engine_type  rand_engine( seed_param.m_seed );

    simulate_engine_type sim_engine( &rand_engine, conf_block );

    config.put_child( CONFIG_BLOCK_K, conf_block );
    writeLog( config, out_path, ".config");

    if( !print_config_only ) {
        unsigned int log_period = (( gen_param.m_size < log_param.m_period ) ? gen_param.m_size : log_param.m_period);

        boost::property_tree::ptree sim_times, stat_times, sim_times_start, sim_times_end;

        timer_type rep_time;
        unsigned int T_gen = gen_param.m_size;
        while( T_gen-- ) {
            timer_type sim_time;
            sim_engine.simulate( );
            sim_time.stop();

            clotho::utility::add_value_array( sim_times, sim_time );
            clotho::utility::add_value_array( sim_times_start, sim_time.getStart() );
            clotho::utility::add_value_array( sim_times_end, sim_time.getStop() );

            if( !( --log_period ) ) {
                boost::property_tree::ptree stat_log;
                timer_type stat_time;
                random_sample_type ss( &rand_engine, sim_engine.getChildPopulation()->haploid_genome_count(), sim_engine.getChildPopulation()->haploid_genome_count() );
                population_analyzer an( &ss, sim_engine.getChildPopulation(), sim_engine.getAlleleSpace() );
                an( stat_log );

                size_t samp_idx = 0;
                BOOST_FOREACH( auto& v, log_param.m_sampling ) {
                    if( v.sample_size <= 0 ) {
                        ++samp_idx;
                        continue;
                    }

                    std::ostringstream oss;
                    oss << "sample." << samp_idx++;

                    random_sample_type samp( &rand_engine, sim_engine.getChildPopulation()->haploid_genome_count(), v.sample_size );

                    population_analyzer _an( &samp, sim_engine.getChildPopulation(), sim_engine.getAlleleSpace() );

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
        perform_log.put_child( "performance.generations.start", sim_times_start );
        perform_log.put_child( "performance.generations.stop", sim_times_end );
        perform_log.put_child( "performance.stats", stat_times );

        sim_engine.getPerformanceResults( perform_log );

        writeLog( perform_log, out_path, ".performance" );
    }

    return 0;
}

population_analyzer::population_analyzer( random_sample_type * s, sequence_space_type * ss, allele_space_type * alls ) :
    m_seqs(ss)
    , m_alleles( alls )
    , m_samp( s )
{}

void population_analyzer::operator()( boost::property_tree::ptree & log ) {
    allele_freq_type af;
//    pairwise_diff_type pd;

    if( m_samp ) {
        af.evaluate( *m_seqs, m_samp->begin(), m_samp->end() );
//        pd.evaluate( *m_samp->getPopulation(), m_samp->begin(), m_samp->end() );
    } else {
        assert( false );
    }

    boost::property_tree::ptree af_log, pd_log;
    af.recordResults( af_log );
//    pd.recordResults( pd_log );

    log.add_child( "allele_frequency", af_log );
//    log.add_child( "pairwise_difference", pd_log );
}
