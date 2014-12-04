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
#include <boost/foreach.hpp>

#include "clotho/powerset/variable_subset_iterator.hpp"
#include "clotho/genetics/individual_phenotyper.hpp"

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

typedef individual_phenotyper< individual_type, no_type >    individual_phenotyper_type;

const string BASE_SEQUENCE_BIAS_K = "base_bias";
const string TRAIT_BLOCK_K = "traits";
const string ALLELE_BLOCK_K = "allele";
const string NEUTRAL_P_K = "neutral.p";

void statsPopulation( population_type * p, allele_set_type * alleles, boost::property_tree::ptree & _log );

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
                statsPopulation(sim.getParentPopulation(), sim.getAlleles(), stat_log );

                individual_phenotyper_type ind_pheno;
                std::vector< typename individual_phenotyper_type::result_type > pop_phenotypes;

                std::transform( sim.getParentPopulation()->begin(), sim.getParentPopulation()->end(), std::back_inserter( pop_phenotypes ), ind_pheno );

                log_type pheno_log;
                BOOST_FOREACH( auto& v, pop_phenotypes ) {
                    clotho::utility::add_value_array( pheno_log, v );
                }

                stat_log.add_child( "phenotypes", pheno_log);

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

void buildPopulationRefMap( population_type * p, ref_map_type & m, allele_dist_type & a ) {
    population_iterator pit = p->begin();
    while( pit != p->end() ) {
        ref_map_iterator rit = m.find( pit->first );
        if( rit == m.end() ) {
            m.insert( std::make_pair( pit->first, 1 ));
        } else {
            ++(rit->second);
        }

        rit = m.find( pit->second);
        if( rit == m.end() ) {
            m.insert( std::make_pair( pit->second, 1 ) );
        } else {
            ++(rit->second);
        }
        ++pit;
    }

    for( ref_map_iterator rit = m.begin(); rit != m.end(); ++rit ) {
        unsigned int n = rit->second;
        size_t idx = rit->first->find_first();
        while( idx !=  sequence_type::bitset_type::npos ) {
            a[ idx ] += n;
            idx = rit->first->find_next( idx );
        }
    }
}

void statsPopulation( population_type * p, allele_set_type * alleles, boost::property_tree::ptree & _log ) {

    typedef std::map< unsigned int, unsigned int >  count_map;
    typedef count_map::iterator                     count_iterator;

    count_map allele_counts;

    ref_map_type seq_ref_counts;
    allele_dist_type allele_dist( alleles->variable_allocated_size(), 0);
    buildPopulationRefMap( p, seq_ref_counts, allele_dist );

    unsigned int nExpSeq = 0;
    accum::accumulator_set< double, accum::stats< accum::tag::weighted_mean, accum::tag::weighted_median, accum::tag::weighted_variance >, double > acc;
    BOOST_FOREACH( auto& v, seq_ref_counts ) {
        double acnt = (( v.first ) ? v.first->count() : 0 );
        double e_rcnt = v.second;

        count_iterator a_it = allele_counts.find( acnt );
        if( a_it == allele_counts.end() ) {
            allele_counts.insert( std::make_pair( acnt, e_rcnt ) );
        } else {
            a_it->second += e_rcnt;
        }

        acc( acnt, accum::weight = e_rcnt );
        nExpSeq += e_rcnt;
    }

    double ave_alleles_per_sequence = accum::weighted_mean(acc);
    double med_allele = accum::weighted_median( acc );
    double var_allele = accum::weighted_variance( acc );

    // for validation purposes
    allele_set_type::cfamily_iterator f_it = alleles->family_begin(), f_end = alleles->family_end();
    unsigned int nNotInPop = 0;
    while( f_it != f_end ) {
        if( seq_ref_counts.find( *f_it++ ) == seq_ref_counts.end()) {
            ++nNotInPop;
        }
    }

    assert( nNotInPop == 0 );

    boost::property_tree::ptree all_dist, all_freq;
    for( count_iterator c_it = allele_counts.begin(); c_it != allele_counts.end(); ++c_it) {
        clotho::utility::add_value_array( all_dist, *c_it );
    }

    typename allele_set_type::cvariable_iterator e_it = alleles->variable_begin();
    for( allele_dist_type::iterator a_it = allele_dist.begin(); a_it != allele_dist.end(); ++a_it, ++e_it ) {
        if( (*a_it) > 0 ) {
            std::ostringstream oss;
            oss << *e_it;
            clotho::utility::add_value_array( all_freq, std::make_pair(oss.str(), (*a_it)));
        } else {
            clotho::utility::add_value_array( all_freq, "");
        }
    }

    _log.put( "population.family_size", alleles->family_size() );

    _log.put( "population.reference_sequences", ((allele_counts.find(0) != allele_counts.end())? allele_counts[0] : 0 ) );

    _log.put( "population.total_sequences", nExpSeq );
    _log.put( "sequences.allele_count_distribution.format", "[alleles per sequence, number of sequences in population]" );
    _log.add_child( "sequences.allele_count_distribution.value", all_dist );
    _log.put( "sequences.alleles_per.mean", ave_alleles_per_sequence);
    _log.put( "sequences.alleles_per.min", allele_counts.begin()->first );
    _log.put( "sequences.alleles_per.max", allele_counts.rbegin()->first );
    _log.put( "sequences.alleles_per.median", med_allele );
    _log.put( "sequences.alleles_per.variance", var_allele );

    _log.put( "alleles.variable_count", alleles->variable_size() );
    _log.put( "alleles.fixed_count", alleles->fixed_size() );
    _log.put( "alleles.free_size", alleles->free_size() );

    _log.put( "alleles.frequencies.format", "[allele, count in population]" );
    _log.add_child( "alleles.frequencies.value", all_freq );
}