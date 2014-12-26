#ifndef SIMULATE_ENGINE_HPP_
#define SIMULATE_ENGINE_HPP_

#include <algorithm>
#include <vector>
#include <map>
#include <unordered_map>
#include <boost/property_tree/ptree.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include "simulation_config.h"

#include "simulate_engine_base.hpp"

#include "clotho/genetics/sequence_mutator.hpp"
#include "clotho/genetics/sequence_generator.hpp"

#include "clotho/genetics/individual_initializer.hpp"
#include "clotho/genetics/individual_selector.hpp"
#include "clotho/genetics/individual_reproduction.hpp"
#include "clotho/genetics/individual_generator.hpp"
#include "clotho/genetics/individual_fitness.hpp"
#include "clotho/genetics/individual_resetter.hpp"
//#include "clotho/genetics/individual_phenotyper.hpp"

#include "clotho/genetics/population_phenotyper.hpp"

#include "infinite_site.hpp"

#include "clotho/genetics/recombiner.hpp"

#include "clotho/powerset/powerset_no_dup_pred.hpp"
#include "clotho/powerset/variable_subset_recombination.hpp"
#include "clotho/powerset/variable_subset_fitness.hpp"

#include "clotho/utility/parameter_space.hpp"

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/weighted_mean.hpp>
#include <boost/accumulators/statistics/weighted_median.hpp>
#include <boost/accumulators/statistics/weighted_variance.hpp>

#include "fitness_generator.hpp"

namespace accum=boost::accumulators;

extern const string SAMPLING_K;

template < class URNG, class AlleleType, class LogType, class TimerType >
class simulate_engine {
public:
    typedef simulate_engine_base< URNG, AlleleType, LogType, TimerType > base_type;

    typedef typename base_type::log_type            log_type;

    typedef typename base_type::rng_type            rng_type;
    typedef typename base_type::allele_type         allele_type;
    typedef typename base_type::allele_set_type     allele_set_type;
    typedef typename base_type::sequence_type       sequence_type;
    typedef typename base_type::sequence_pointer    sequence_pointer;
    typedef typename base_type::individual_type     individual_type;

    typedef typename base_type::allele_generator    allele_generator;

    // recombination typedefs
    typedef typename base_type::classifier_type     classifier_type;
    typedef clotho::recombine::recombination< sequence_type, classifier_type >  recombination_engine_type;
    typedef recombiner< sequence_type, recombination_engine_type >              recombination_method;
    typedef clotho::utility::random_generator< rng_type, recombination_method > recombination_method_generator;
    typedef typename base_type::population_type                                 population_type;
    typedef typename population_type::iterator                                  population_iterator;

    typedef clotho::utility::random_generator< rng_type, infinite_site< sequence_type > >       mutation_generator_type;
    typedef sequence_generator< sequence_pointer >                                              sequence_generator_type;
    typedef sequence_mutator< sequence_type, mutation_generator_type >                          sequence_mutator_type;
    typedef clotho::utility::random_generator< rng_type, sequence_mutator_type >                sequence_mutator_generator;
    typedef individual_initializer< individual_type, sequence_generator_type >   individual_initializer_type;
    typedef individual_selector< rng_type >             individual_selector_type;
    typedef individual_reproduction< individual_type
                , sequence_mutator_generator
                , recombination_method_generator > individual_reproduction_type;
    typedef individual_generator< population_type, individual_selector_type, individual_reproduction_type >     individual_generator_type;

    typedef individual_resetter< individual_type >                  individual_resetter_type;
//    typedef individual_phenotyper< individual_type, no_type >       individual_phenotyper_type;
//
//    typedef std::vector< typename individual_phenotyper_type::result_type > population_phenotypes;

    typedef population_phenotyper< population_type >                population_phenotyper_type;
    typedef typename population_phenotyper_type::phenotype_type     phenotype_type;
    typedef typename population_phenotyper_type::result_type        population_phenotypes;

    // fitness typedefs
//    typedef clotho::fitness::fitness_method< double, clotho::fitness::multiplicative_heterozygous_tag > het_fit_type;
//    typedef clotho::fitness::fitness_method< double, clotho::fitness::multiplicative_homozygous_tag >    alt_hom_fit_type;
//    typedef clotho::fitness::no_fit                                             ref_hom_fit_type;
//    typedef clotho::fitness::fitness< sequence_type, het_fit_type, alt_hom_fit_type, ref_hom_fit_type, double > fitness_type;
//    typedef typename fitness_type::result_type                                  fitness_result_type;
//    typedef individual_fitness< fitness_type >                                  fitness_operator;
//
#ifdef CONSTFIT
    typedef constant_fitness_generator  fitness_generator_type;
#elif defined QFIT
    typedef quadratic_fitness_generator fitness_generator_type;
#else
    typedef normal_fitness_generator    fitness_generator_type;
#endif  // QFIT


    typedef typename fitness_generator_type::result_type                        fitness_operator;
    typedef typename fitness_operator::result_type                              fitness_result_type;
    typedef std::vector< fitness_result_type >                                  population_fitnesses;

    // statistic typedefs
    typedef std::map< sequence_pointer, unsigned int >  ref_map_type;
    typedef typename ref_map_type::iterator             ref_map_iterator;
    typedef std::vector< unsigned int >                 allele_dist_type;

    simulate_engine( boost::property_tree::ptree & config ) :
        m_rng()
        , m_founder_size( DEFAULT_POPULATION_SIZE )
        , m_seq_mut_gen( m_rng, config )
        , m_rec_met_gen( m_rng, config )
        , m_repro( m_seq_mut_gen, m_rec_met_gen )
        , m_fit_gen( config )
        , m_parent( &m_pop )
        , m_child( &m_buffer )
        , m_parent_pheno( &m_pheno_buff1 )
        , m_child_pheno( &m_pheno_buff2 )
        , m_parent_fit( &m_fit_buff1 )
        , m_child_fit( &m_fit_buff2 )
    {
        parseConfig( config );
        initialize();
    }

    void simulate() {
        simulate( m_founder_size );
    }

    void simulate( unsigned int p_size ) {
/*
        population_fitnesses pfit;
        pfit.reserve( m_parent->size() );

        fitness_type fit;
        fitness_operator fit_op( m_alleles, fit );

//        fitness_result_type tot_fit = 0., exp_fit = 0.;
        if( !fit_op.selected_count() ) {
            std::fill_n( std::back_inserter(pfit), m_parent->size(), 1. );
//            tot_fit = (fitness_result_type)m_parent->size();
//            exp_fit = 1.;
        } else {
            std::transform( m_parent->begin(), m_parent->end(), std::back_inserter(pfit), fit_op );
//            tot_fit = fit_op.total_fitness();
//            exp_fit = fit_op.expected_fitness();
        }
*/
//       individual_selector_type sel( m_rng, pfit.begin(), pfit.end() );

        individual_selector_type sel( m_rng, m_parent_fit->begin(), m_parent_fit->end() );
        individual_generator_type ind_gen( m_parent, sel, m_repro );

        std::generate_n( std::back_inserter( *m_child ), p_size, ind_gen );

        updatePhenotypes( *m_child_pheno, *m_child );
        updateFitness( *m_child_fit, *m_child_pheno );
    }

    void reset_parent() {
        m_parent->clear();

        std::swap( m_parent, m_child );
        std::swap( m_parent_pheno, m_child_pheno );
        std::swap( m_parent_fit, m_child_fit );

        m_alleles.pruneSpace();

        m_child_pheno->clear();
        m_child_fit->clear();
    }

    population_type *   getParentPopulation() { return m_parent; }
    population_type *   getChildPopulation() { return m_child; }

    allele_set_type *   getAlleles() { return &m_alleles; }

    log_type &          getLog() { return m_log; }

    void                clearLog() { m_log.clear(); }

    log_type            getState() {
        log_type state;

        log_type p, c, a;
//        state_of< population_type >::record( *m_parent, p );
//        state_of< population_type >::record( *m_child, c );
//        state_of< allele_set_type >::record( m_alleles, a );

        state.put( "population.parent", p );
        state.put( "population.child", c );

        state.put( "alleles", a);
        return state;
    }

    void computeStats( boost::property_tree::ptree & log ) {
        boost::property_tree::ptree pheno_log;
        BOOST_FOREACH( auto& v, *m_parent_pheno ) {
            clotho::utility::add_value_array( pheno_log, v );
        }

        log.add_child( "phenotypes", pheno_log);

        boost::property_tree::ptree fit_log;
        BOOST_FOREACH( auto& v, *m_parent_fit ) {
            clotho::utility::add_value_array( fit_log, v );
        }

        log.add_child( "fitness", fit_log );

        population_statistics( m_parent, &m_alleles, log );
        unsigned int i = 0;
        BOOST_FOREACH( auto& v, m_sampling_sizes ) {
            boost::property_tree::ptree l;
            random_sampler( m_parent, m_parent_pheno, m_parent_fit, v, l );

            l.put("size", v );

            std::ostringstream oss;
            oss << "sample." << (i++);

            log.add_child( oss.str(), l );
        }
    }

protected:

    inline void updatePhenotypes( population_phenotypes & phenos, population_type & p ) {
        population_phenotyper_type ppheno;
        BOOST_FOREACH( auto& i, p ) {
            phenos.push_back( ppheno(i) );
        }
    }

    inline void updateFitness( population_fitnesses & fit, population_phenotypes & phenos ) {
        fitness_operator pfit = m_fit_gen( phenos.begin(), phenos.end() );
        BOOST_FOREACH( auto& p, phenos ) {
            fit.push_back( pfit(p) );
        }
    }

    void parseConfig( boost::property_tree::ptree & config ) {
        std::ostringstream oss;
        oss << CONFIG_BLOCK_K << "." << RNG_BLOCK_K << "." << SEED_K;

        if( config.get_child_optional( oss.str() ) == boost::none ) {
            config.put( oss.str(), 0 );
        } else {
            unsigned int seed = config.get< unsigned int >( oss.str(), 0 );
            m_rng.seed( seed );
        }

        oss.str("");
        oss.clear();

        oss << CONFIG_BLOCK_K << "." << POP_BLOCK_K << "." << SIZE_K;
        if( config.get_child_optional( oss.str() ) == boost::none ) {
            config.put( oss.str(), m_founder_size );
        } else {
            m_founder_size = config.get< unsigned int >(oss.str(), m_founder_size );
        }

        oss.str("");
        oss.clear();

        oss << CONFIG_BLOCK_K << "." << LOG_BLOCK_K << "." << SAMPLING_K;

        if( config.get_child_optional( oss.str() ) == boost::none ) {
            config.put( oss.str(), "" );
        } else {
            BOOST_FOREACH( auto& v, config.get_child( oss.str() ) ) {
                std::ostringstream tmp;
                tmp << v.second.data();

                unsigned int s = boost::lexical_cast< unsigned int >( tmp.str() );
                m_sampling_sizes.push_back(s);
            }
        }
    }

    void initialize( ) {
        m_pop.clear();
        m_buffer.clear();

        m_pop.reserve( m_founder_size );
        m_buffer.reserve( m_founder_size );

        // results in the existence of multiple 'empty' sequences in the family
        sequence_generator_type sgen( m_alleles );
        individual_initializer_type igen( sgen );
        std::generate_n( std::back_inserter( m_pop ), m_founder_size, igen );

        //sequence_pointer empty_seq = m_alleles.create_subset();

        //std::fill_n( std::back_inserter( *m_parent ), m_founder_size, empty_seq );
        std::fill_n( std::back_inserter( *m_parent_pheno ), m_founder_size, phenotype_type(1, 0.) );
        updateFitness( *m_parent_fit, *m_parent_pheno );
    }

    void random_sampler( population_type * p, population_phenotypes * phenos, population_fitnesses * fits, unsigned int size, log_type & l ) {
        population_type sub_pop;
        boost::property_tree::ptree pheno_log, fit_log;
        boost::random::uniform_int_distribution< unsigned int > ind_rand(0, p->size() - 1 );
        for( unsigned int i = 0; i < size; ++i ) {
            unsigned int idx = ind_rand( m_rng );
            sub_pop.push_back( p->at(idx) );
            clotho::utility::add_value_array( pheno_log, phenos->at(idx) );
            clotho::utility::add_value_array( fit_log, fits->at(idx) );
        }

        l.add_child( "phenotypes", pheno_log);
        l.add_child( "fitness", fit_log );

        population_statistics( &sub_pop, &m_alleles, l );
    }

    void buildRefMap( population_type * p, ref_map_type & m, allele_dist_type & a ) {

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

    void population_statistics( population_type * p, allele_set_type * alleles, boost::property_tree::ptree & _log ) {

        typedef std::map< unsigned int, unsigned int >  count_map;
        typedef count_map::iterator                     count_iterator;

        count_map allele_counts;

        ref_map_type seq_ref_counts;
        allele_dist_type allele_dist( alleles->variable_allocated_size(), 0);
        buildRefMap(p, seq_ref_counts, allele_dist );

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

        //_log.put( "population.family_size", seq_ref_counts.size() );
        //_log.put( "population.reference_sequences", ((allele_counts.find(0) != allele_counts.end())? allele_counts[0] : 0 ) );
        //_log.put( "population.total_sequences", nExpSeq );
        //
        _log.put( "population.sequences.unique", seq_ref_counts.size() );
        _log.put( "population.sequences.num_empty", ((allele_counts.find(0) != allele_counts.end() )?allele_counts[0] : 0));
        _log.put( "population.size", p->size() );

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

    rng_type        m_rng;
    unsigned int    m_founder_size;
    sequence_mutator_generator m_seq_mut_gen;
    recombination_method_generator m_rec_met_gen;

    individual_reproduction_type    m_repro;
    fitness_generator_type          m_fit_gen;

    log_type            m_log;
    population_type     m_pop, m_buffer;
    population_type     * m_parent, * m_child;

    population_phenotypes   m_pheno_buff1, m_pheno_buff2;
    population_phenotypes   * m_parent_pheno, * m_child_pheno;

    population_fitnesses m_fit_buff1, m_fit_buff2;
    population_fitnesses * m_parent_fit, * m_child_fit;

    allele_set_type m_alleles;

    std::vector< unsigned int > m_sampling_sizes;
};

namespace clotho {
namespace utility {

template < class URNG, class AlleleType, class LogType, class TimerType >
struct parameter_space< simulate_engine< URNG, AlleleType, LogType, TimerType > > {

    static void build_parameters( boost::property_tree::ptree & params ) {
        std::ostringstream oss;
        oss << CONFIG_BLOCK_K << "." << POP_BLOCK_K << "." << SIZE_K;
        params.put( oss.str() + ".value", 1000 );
        params.put( oss.str() + ".type", "uint");
        params.put( oss.str() + ".description", "Founder Population Size" );
    }
};

}   // namespace utility
}   // namespace clotho

#endif  // SIMULATE_ENGINE_HPP_
