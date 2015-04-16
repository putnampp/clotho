#include "config.h"
#include "qtl_config.hpp"

#include <iostream>
#include <algorithm>
#include <map>

#include "clotho/utility/log_helper.hpp"
#include "common_commandline.h"

#include "clotho/genetics/qtl_allele.h"
#include "clotho/utility/timer.hpp"
#include "simulate_engine.hpp"

#include <boost/random/mersenne_twister.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
//#include <boost/foreach.hpp>
//
//#include "clotho/powerset/variable_subset_iterator.hpp"
//#include "clotho/genetics/individual_phenotyper.hpp"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/weighted_mean.hpp>
#include <boost/accumulators/statistics/weighted_median.hpp>
#include <boost/accumulators/statistics/weighted_variance.hpp>

#include "clotho/genetics/fitness_toolkit.hpp"

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

const string BASE_SEQUENCE_BIAS_K = "base_bias";
const string TRAIT_BLOCK_K = "traits";
const string ALLELE_BLOCK_K = "allele";
const string NEUTRAL_P_K = "neutral.p";

const string SAMPLING_K = "sampling_size";

//const string FITNESS_BLOCK_K = "fitness";
const string QUADRATIC_SCALE_K = "quadratic.scale";
const string CONSTANT_K = "constant";

const string SIZE_K = "size";
const string PAIRWISE_K = "pairwise";

const string ENGINE_BLOCK_K = "engine";
const string POWERSET_BLOCK_K = "powerset";
const string SUBSET_BLOCK_K = "subset";
const string CLASSIFIER_BLOCK_K = "classifier";
const string TYPE_K = "type";

void get_engine_config( const std::string & out_path );

int main( int argc, char ** argv ) {

    log_type config;
    int res = parse_commandline(argc, argv, config);
    if( res ) {
        return res;
    }

    log_type conf_child = ( (config.get_child_optional( CONFIG_BLOCK_K) == boost::none) ? config : config.get_child( CONFIG_BLOCK_K ));

    const unsigned int nRep = conf_child.get< unsigned int >( REPEAT_K, 1 );
    const unsigned int nGen = conf_child.get< unsigned int >( GEN_BLOCK_K + "." + SIZE_K, 1);
    const unsigned int nLog = conf_child.get< unsigned int >( LOG_BLOCK_K + "." + PERIOD_K, -1);
    const unsigned int seed = conf_child.get< unsigned int >( RNG_BLOCK_K + "." + SEED_K, 0 );
    string out_path = conf_child.get<string>( OUTPUT_K, "");

    get_engine_config( out_path );

    rng_type rng(seed);

    for( unsigned int i = 0; i < nRep; ++i ) {
        // change the seed value of the random number generator
        rng.discard( 15 );
        const unsigned int tmp_seed = rng();
        log_type rep_child_conf = conf_child;

        rep_child_conf.put( RNG_BLOCK_K + "." + SEED_K, tmp_seed );

        simulate_type sim( rep_child_conf );

        if( nGen == 0 ) {
            fitness_toolkit::getInstance()->tool_configurations( rep_child_conf );

            log_type tmp;
            tmp.add_child( CONFIG_BLOCK_K, rep_child_conf );
            boost::property_tree::write_json( std::cerr, tmp );

            continue;
        }

        unsigned int log_period = ((nGen < nLog) ? nGen : nLog);
        log_type sim_times, stat_times;

        timer_type rep_time;
        for( unsigned int j = 0; j < nGen; ++j ) {
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

                log_period = ((j + nLog < nGen) ? nLog : (nGen - j - 1) );

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

void get_engine_config( const std::string & out_path ) {
    boost::property_tree::ptree compile_log;

    compile_log.put( ENGINE_BLOCK_K + ".description", "Simulator compiled objects; READ ONLY");

    compile_log.put(ENGINE_BLOCK_K + "." + POWERSET_BLOCK_K + "." + SUBSET_BLOCK_K  +"." + TYPE_K, STR( SUBSETTYPE ) );
    compile_log.put(ENGINE_BLOCK_K + "." + REC_BLOCK_K + "." + TYPE_K, STR( RECOMBTYPE ) );
    compile_log.put(ENGINE_BLOCK_K + "." + POWERSET_BLOCK_K +"." + SIZE_K, BLOCK_UNIT_SIZE );

    compile_log.put( ENGINE_BLOCK_K + "." + POWERSET_BLOCK_K + "." + SUBSET_BLOCK_K + "." + REC_BLOCK_K + ".tag0", STR( RECOMBINE_INSPECT_METHOD ) );
    compile_log.put( ENGINE_BLOCK_K + "." + POWERSET_BLOCK_K + "." + SUBSET_BLOCK_K + "." + REC_BLOCK_K + ".tag1", STR( BIT_WALK_METHOD ) );

    compile_log.put( ENGINE_BLOCK_K + ".reproduction_method.type", STR(REPRODUCTION_METHOD_TAG));

    if( out_path.empty() ) {
        boost::property_tree::write_json( std::cout, compile_log );
    } else {
        std::ostringstream oss;
        oss << out_path << ".engine_compilation.json";

        boost::property_tree::write_json( oss.str(), compile_log );
    }
}
