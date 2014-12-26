#include "config.h"

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

//typedef individual_phenotyper< individual_type, no_type >    individual_phenotyper_type;

const string BASE_SEQUENCE_BIAS_K = "base_bias";
const string TRAIT_BLOCK_K = "traits";
const string ALLELE_BLOCK_K = "allele";
const string NEUTRAL_P_K = "neutral.p";

const string SAMPLING_K = "sampling_size";

const string FITNESS_BLOCK_K = "fitness";
const string QUADRATIC_SCALE_K = "quadratic.scale";
const string CONSTANT_K = "constant";
/*
void random_sample( rng_type * rng, population_type * p, allele_set_type * alleles, unsigned int nSamples, boost::property_tree::ptree & _log );
void statsPopulation( population_type * p, allele_set_type * alleles, boost::property_tree::ptree & _log );
*/
int main( int argc, char ** argv ) {

    log_type config;
    int res = parse_commandline(argc, argv, config);
    if( res ) {
        return res;
    }
    
    const unsigned int nRep = config.get< unsigned int >( CONFIG_BLOCK_K + "." + REPEAT_K, 1 );
    const unsigned int nGen = config.get< unsigned int >( CONFIG_BLOCK_K + "." + GEN_BLOCK_K + "." + SIZE_K, 1);
    const unsigned int nLog = config.get< unsigned int >( CONFIG_BLOCK_K + "." + LOG_BLOCK_K + "." + PERIOD_K, -1);
    const unsigned int seed = config.get< unsigned int >( CONFIG_BLOCK_K + "." + RNG_BLOCK_K + "." + SEED_K, 0 );

    rng_type rng(seed);

    string out_path = config.get<string>( CONFIG_BLOCK_K + "." + OUTPUT_K, "");

    for( unsigned int i = 0; i < nRep; ++i ) {
        // change the seed value of the random number generator
        rng.discard( 15 );
        const unsigned int tmp_seed = rng();
        config.put( CONFIG_BLOCK_K + "." + RNG_BLOCK_K + "." + SEED_K, tmp_seed );

        simulate_type sim( config );

        unsigned int log_period = ((nGen < nLog) ? nGen : nLog);
        for( unsigned int j = 0; j < nGen; ++j ) {
            sim.simulate();

            sim.reset_parent();
            if( !(--log_period) ) {
                log_type stat_log;
                sim.computeStats( stat_log );

                // combine simulation log and configuration log into single object
                BOOST_FOREACH( auto& upd, sim.getLog() ) {
                    stat_log.put_child( upd.first, upd.second );
                }
                sim.clearLog();

                BOOST_FOREACH( auto& upd, config ) {
                    stat_log.put_child(upd.first, upd.second );
                }

                log_period = ((j + log_period < nGen) ? nLog : (nGen - j - 1) );

                if( !stat_log.empty() ) {
                    if( out_path.empty() ) {
                        boost::property_tree::write_json( std::cout, stat_log );
                    } else {
                        std::ostringstream oss;
                        oss << out_path << "." << i << "." << j << ".json";

                        boost::property_tree::write_json( oss.str(), stat_log );
                    }
                }
            }
        }
    }

    return 0;
}
