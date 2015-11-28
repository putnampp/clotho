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
#include "config.h"
#include "qtl_config.hpp"

#include <iostream>
#include <algorithm>
#include <map>

#include "clotho/utility/clotho_strings.hpp"
#include "clotho/utility/log_helper.hpp"
#include "common_commandline.h"

#include "clotho/genetics/qtl_allele.h"
#include "clotho/utility/timer.hpp"
#include "simulate_engine.hpp"

#include <boost/random/mersenne_twister.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/weighted_mean.hpp>
#include <boost/accumulators/statistics/weighted_median.hpp>
#include <boost/accumulators/statistics/weighted_variance.hpp>

#include "clotho/fitness/fitness_toolkit.hpp"

#include "../generation_parameter.hpp"
#include "clotho/random/seed_parameter.hpp"
#include "qtl_logging_parameter.hpp"

#define XSTR( x ) #x
#define STR( x )  XSTR( x )

namespace accum=boost::accumulators;

typedef boost::random::mt19937                                              rng_type;
typedef qtl_allele                                                          allele_type;
typedef boost::property_tree::ptree                                         log_type;
typedef clotho::utility::timer                                              timer_type;

typedef simulate_engine< rng_type, allele_type, log_type, timer_type > simulate_type;

// allele set typedef
typedef typename simulate_type::sequence_type    sequence_type;
typedef typename simulate_type::sequence_pointer sequence_pointer;
typedef typename simulate_type::allele_set_type allele_set_type;

// population typedefs
typedef typename simulate_type::individual_type                        individual_type;
typedef typename simulate_type::population_type                        population_type;
typedef typename population_type::iterator                             population_iterator;

typedef std::map< sequence_pointer, unsigned int >  ref_map_type;
typedef typename ref_map_type::iterator             ref_map_iterator;
typedef std::vector< unsigned int >                 allele_dist_type;

//const string BASE_SEQUENCE_BIAS_K = "base_bias";
//const string TRAIT_BLOCK_K = "traits";
//const string ALLELE_BLOCK_K = "allele";
//const string NEUTRAL_P_K = "neutral.p";

//const string SAMPLING_K = "sampling_size";

//const string CONSTANT_K = "constant";

//const string SIZE_K = "size";
//const string PAIRWISE_K = "pairwise";

const string ENGINE_BLOCK_K = "engine";
const string POWERSET_BLOCK_K = "powerset";
const string SUBSET_BLOCK_K = "subset";
const string CLASSIFIER_BLOCK_K = "classifier";
//const string TYPE_K = "type";

void write_engine_config( const std::string & out_path );

int main( int argc, char ** argv ) {

    log_type config;
    int res = parse_commandline(argc, argv, config);
    if( res ) {
        return res;
    }

    log_type conf_child;
    conf_child = config.get_child( CONFIG_BLOCK_K, conf_child );

    const unsigned int nRep = conf_child.get< unsigned int >( REPETITION_K, 1 );

    generation_parameter gen_param( conf_child );

    qtl_logging_parameter log_param( conf_child );
    seed_parameter<  > seed_param( conf_child );

    std::cout << "Logging period: " << log_param.m_period << std::endl;

    string out_path = conf_child.get<string>( PREFIX_K, "" );

    write_engine_config( out_path );

    rng_type rng(seed_param.m_seed);

    for( unsigned int i = 0; i < nRep; ++i ) {
        // change the seed value of the random number generator
        rng.discard( 15 );
        const unsigned int tmp_seed = rng();
        log_type rep_child_conf = conf_child;

        rep_child_conf.put( RNG_BLOCK_K + "." + SEED_K, tmp_seed );

        simulate_type sim( rep_child_conf );

        if( gen_param.m_size < 1 ) {
            fitness_toolkit::getInstance()->tool_configurations( rep_child_conf );
            population_growth_toolkit::getInstance()->tool_configurations( rep_child_conf );

            log_type tmp;
            tmp.add_child( CONFIG_BLOCK_K, rep_child_conf );
            if( out_path.empty() ) {
                boost::property_tree::write_json( std::cerr, tmp );
            } else {
                std::ostringstream oss;
                oss << out_path << "_config.json";
                boost::property_tree::write_json( oss.str(), tmp );
            }

            continue;
        }

        unsigned int log_period = ((gen_param.m_size < log_param.m_period) ? gen_param.m_size : log_param.m_period);
        log_type sim_times, stat_times;

        timer_type rep_time;
        for( unsigned int j = 0; j < gen_param.m_size; ++j ) {
            timer_type sim_time;
            sim.simulate(j);
            sim_time.stop();

            clotho::utility::add_value_array( sim_times, sim_time );

            sim.reset_parent();
            if( !(--log_period) ) {
                timer_type stat_time;

                log_type stat_log;
                sim.computeStats( stat_log );

                // combine simulation log and configuration log into single object
                BOOST_FOREACH( auto& upd, sim.getLog() ) {
                    stat_log.put_child( upd.first, upd.second );
                }
                sim.clearLog();

//                BOOST_FOREACH( auto& upd, config ) {
//                    stat_log.put_child(upd.first, upd.second );
//                }
                stat_log.put_child( CONFIG_BLOCK_K, rep_child_conf);

                log_period = ((j + log_param.m_period < gen_param.m_size) ? log_param.m_period : (gen_param.m_size - j - 1) );

                if( !stat_log.empty() ) {
                    if( out_path.empty() ) {
                        boost::property_tree::write_json( std::cout, stat_log );
                    } else {
                        std::ostringstream oss;
                        oss << out_path << "." << i << "." << j << ".json";

                        boost::property_tree::write_json( oss.str(), stat_log );
                    }
                }
                stat_time.stop();
                clotho::utility::add_value_array( stat_times, stat_time );
            }
        }
        rep_time.stop();

        log_type perform_log;
        perform_log.add( "performance.runtime", rep_time.elapsed().count() );
        perform_log.put_child( "performance.simulate", sim_times );
        perform_log.put_child( "performance.stats", stat_times );

        if( out_path.empty() ) {
            boost::property_tree::write_json( std::cout, perform_log );
        } else {
            std::ostringstream oss;
            oss << out_path << "." << i << ".performance.json";

            boost::property_tree::write_json( oss.str(), perform_log );
        }
    }

    return 0;
}

void write_engine_config( const std::string & out_path ) {

    boost::property_tree::ptree sset, recomb, pset;
    recomb.put( "tag0", STR( RECOMBINE_INSPECT_METHOD ) );
    recomb.put( "tag1", STR( BIT_WALK_METHOD ) );
    sset.put_child( "recombine", recomb );
    sset.put( TYPE_K, STR( SUBSETTYPE ) );

    pset.put_child( SUBSET_BLOCK_K, sset );
    pset.put( "block_size", BLOCK_UNIT_SIZE );
//    compile_log.put(ENGINE_BLOCK_K + "." + REC_BLOCK_K + "." + TYPE_K, STR( (RECOMBTYPE) ) );
//    compile_log.put(ENGINE_BLOCK_K + "." + POWERSET_BLOCK_K +"." + SIZE_K, BLOCK_UNIT_SIZE );

//    compile_log.put(ENGINE_BLOCK_K + "." + POWERSET_BLOCK_K + "." + SUBSET_BLOCK_K  +"." + TYPE_K, STR( SUBSETTYPE ) );
//    compile_log.put( ENGINE_BLOCK_K + "." + POWERSET_BLOCK_K + "." + SUBSET_BLOCK_K + "." + REC_BLOCK_K + ".tag0", STR( RECOMBINE_INSPECT_METHOD ) );
//    compile_log.put( ENGINE_BLOCK_K + "." + POWERSET_BLOCK_K + "." + SUBSET_BLOCK_K + "." + REC_BLOCK_K + ".tag1", STR( BIT_WALK_METHOD ) );

//    compile_log.put( ENGINE_BLOCK_K + ".reproduction_method.type", STR(REPRODUCTION_METHOD_TAG));
//
//    compile_log.put( ENGINE_BLOCK_K + ".individual_selector.type", STR(IND_SELECT) );

    boost::property_tree::ptree eng;
    eng.put_child( POWERSET_BLOCK_K, pset );
    eng.put( REC_BLOCK_K + ".type", STR( (RECOMBTYPE) ) );
    eng.put( "reproduction_method.type", STR(REPRODUCTION_METHOD_TAG));
    eng.put( "individual_selector.type", STR(IND_SELECT) );
    eng.put( "description", "Simulator compiled objects; READ ONLY");

    boost::property_tree::ptree compile_log;
    compile_log.put_child( ENGINE_BLOCK_K, eng );

    if( out_path.empty() ) {
        boost::property_tree::write_json( std::cout, compile_log );
    } else {
        std::ostringstream oss;
        oss << out_path << ".engine_compilation.json";

        boost::property_tree::write_json( oss.str(), compile_log );
    }
}
