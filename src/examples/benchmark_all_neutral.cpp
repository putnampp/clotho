//#define LOGGING 1

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
typedef boost::property_tree::ptree state_log_type;
#ifdef LOGGING
state_log_type global_log;

#include "utility/lowest_bit.h"
#endif

#include "common_commandline.h"
#include <cstdlib>
#include <functional>
#include <algorithm>

#include <vector>
#include <memory>
#include <limits>
#include <cmath>

#include "clotho/utility/timer.hpp"

#include "object/individual.hpp"
#include "object/individual_properties.hpp"

#include "genetics/locus_bitset.h"

#include "rng/rng.hpp"

#define _PLOIDY 2

using std::shared_ptr;
using clotho::utility::timer;

#include "reproduce_with_recombination.hpp"
#include "fitness_methods.hpp"
#include "clotho/selection.hpp"
#include "specializations.tcc"

#include "population_graphs.h"

typedef locus_bitset gamete_type;
typedef typename gamete_type::alphabet_t    alphabet_type;
typedef typename gamete_type::allele_type   allele_type;

typedef allele_type *   allele_pointer;
typedef gamete_type *   gamete_pointer;

typedef reproduction::models::mutation::mutate_site< allele_type, alphabet_type, gamete_type >   mmodel_type;
typedef reproduction::models::recombination::no_recomb< _PLOIDY >     rcmodel_type;
typedef reproduction::IndividualReproduction< mmodel_type, rcmodel_type > rmodel_type;

typedef TIndividual< gamete_type, _PLOIDY, system_id, individual_props< gamete_type, _PLOIDY, system_id >> individual_type;
typedef individual_type  *                  individual_pointer;
typedef std::vector< individual_pointer >   environment_type;

struct simulation_configuration {
    state_log_type log;

    unsigned int                nGen, nPop;

    double                      mu, rho;

    gsl_rng *                   my_rng;

    std::shared_ptr< iRNG >     rng;

    simulation_configuration( ) :
        nGen(0), nPop(0), mu( 0.0 ), rho(0.0), my_rng( NULL )
    {
        init_rng();
    }

    simulation_configuration( unsigned int gen, unsigned int pop, double m, double r ) :
        nGen( gen ), nPop(pop), mu(m), rho(r)
    {
        init_rng();
    }

    void init_rng( unsigned int seed = 0 ) {
        gsl_rng_env_setup();
        const gsl_rng_type * T = gsl_rng_default;
        string m_type = T->name;
        unsigned int m_seed = gsl_rng_default_seed;

        my_rng = gsl_rng_alloc( T );
        gsl_rng_set( my_rng, m_seed );

        rng.reset( new GSL_RNG( my_rng, m_type, m_seed ));
    }
};

int main( int argc, char ** argv ) {

    po::variables_map vm;
    if( parse_commandline( argc, argv, vm ) != 0 ) {
        return 0;
    }

    timer runtime_timer;
    state_log_type  log;

    unsigned int nGen = vm[ GENERATIONS_K ].as< unsigned int >();

    gsl_rng_env_setup();
    const gsl_rng_type * T = gsl_rng_default;
    string m_type = T->name;
    unsigned int m_seed = gsl_rng_default_seed;

    gsl_rng * my_rng = gsl_rng_alloc( T );
    gsl_rng_set( my_rng, m_seed );

    shared_ptr< iRNG > rng( new GSL_RNG( my_rng, m_type, m_seed ));
    RandomProcess::initialize( rng );

    double mu = vm[ MUTATION_RATE_K ].as<double>();
    double rho = vm[ RECOMBINATION_RATE_K ].as< double >();

    log.put( "configuration.random_number.generator", rng->getType() );
    log.put( "configuration.random_number.seed", rng->getSeed() );

    log.put( "configuration.generations", nGen );
    log.put( "configuration.population.size", vm[ FOUNDER_SIZE_K ].as< unsigned int >() );
    log.put( "configuration.mutation.model", "infinite site");
    log.put( "configuration.mutation.rate_per_region", mu );

    log.put( "configuration.recombination.model", "uniform random" );
    log.put( "configuration.recombination.rate_per_region", rho );

    log.put( "configuration.region.mutation_rate_per_base", pow( 10.0, -8.0) );
    log.put( "configuration.region.base_per_region", mu / pow( 10.0, -8.0) );
    log.put( "configuration.region.regions_per_individual",  _PLOIDY );

    mmodel_type::initialize( mu, false);

    environment_type population;
    environment_type buffer;

    system_id blank_id;

    fitness_multiplicative< het_fitness< allele_type >, hom_fitness< allele_type >, gamete_type, individual_pointer > fmult;

    timer init_time;
    for( unsigned int i = 0; i < vm[ FOUNDER_SIZE_K ].as< unsigned int >(); ++i) {
        population.push_back( new individual_type() );

        population.back()->getProperties()->inheritFrom( blank_id,  gamete_type::EMPTY.copy() );
        population.back()->getProperties()->inheritFrom( blank_id,  gamete_type::EMPTY.copy() );

        buffer.push_back( new individual_type() );
    }

    std::cerr << gamete_type::EMPTY.copies() << std::endl;
    assert( gamete_type::EMPTY.copies() - 1 == 2 * vm[ FOUNDER_SIZE_K ].as<unsigned int>() );

    init_time.stop();

    environment_type * parent = &population, * child = &buffer;

    unsigned int fitness_size = 0;
    double * fitness = NULL;

    ReproduceWithRecombination<gamete_type, alphabet_type, individual_pointer> repro( mu, rho );

    state_log_type gen_log;
    state_log_type gen_run_log;
    state_log_type gen_size_log;
    state_log_type gen_alpha_log;
    state_log_type gen_seq_log;
    state_log_type gen_fit_log, gen_repro_log, gen_reset_log;

    size_t nSelfing = 0;
    timer sim_time;

    std::ostringstream oss;
    size_t nBlocks = 0;

    for( SystemClock::vtime_t i = 0; i < nGen; ++i ) {
        timer gen_timer;
        timer fit_timer;
        assert( parent != child );
        if( fitness_size < parent->size() ) {
            if( fitness != NULL ) {
                delete [] fitness;
            }
            fitness_size = parent->size();
            fitness = new double[ fitness_size ];
        }

        memset( fitness, 0, sizeof(double) * fitness_size );

#ifdef LOGGING
        {
            std::string k = oss.str();

            state_log_type tmp;
            alphabet_type::getInstance()->logState( tmp );

            global_log.add_child( oss.str(), tmp );
        }
#endif

        // measure fitness of parent population
        //
        double * tmp = fitness;
//        double e_fitness = 0.0;
        for( environment_type::iterator it = parent->begin(); it != parent->end(); it++ ) {
            (*tmp) = 1.0;
            (*tmp) = fmult( (*tmp), (*it) );
            //e_fitness += (*tmp++);
            ++tmp;
        }
        //e_fitness /= (double)parent->size();
        fit_timer.stop();

        timer repro_timer;
        //
        // mate
        //
        clotho::selection::DiscreteSelection< individual_pointer, environment_type * > ds( my_rng, fitness, parent->size() );
        unsigned int child_idx = 0;
        while( child_idx < child->size()) {
            std::pair< individual_pointer, individual_pointer > mate_pair = ds( parent );

            nSelfing += (( mate_pair.first == mate_pair.second ) ? 1 : 0);
            gamete_pointer g = repro( mate_pair.first, i );
            (*child)[child_idx]->getProperties()->inheritFrom(blank_id, g);

            gamete_pointer g1 = repro( mate_pair.second, i );
            (*child)[child_idx]->getProperties()->inheritFrom(blank_id, g1);

            (*child)[child_idx++]->getProperties()->setDOB( i );
        }

        repro_timer.stop();

        timer reset_timer;
        assert( parent->size() == vm[ FOUNDER_SIZE_K ].as< unsigned int > () );
        for( environment_type::iterator it = parent->begin(); it != parent->end(); it++ ) {
            (*it)->reset();
        }

        std::swap( parent, child );

        nBlocks = locus_bitset::updateActiveAlphabet();
        reset_timer.stop();

        gen_timer.stop();

        size_t parent_pop_seq_count = locus_bitset::activeCount();
        size_t parent_pop_mut_count = alphabet_type::getInstance()->active_count(); // performs dynamic_bitset::count(); linear computation of bitset Hamming Weight
        
        state_log_type a, g, r, s, t, f, m, c;
        a.put("", parent_pop_mut_count );
        gen_alpha_log.push_back( std::make_pair("", a ) );

        g.put("", i );
        gen_log.push_back( std::make_pair("", g ) );

        r.put("", gen_timer.elapsed().count() );
        gen_run_log.push_back( std::make_pair("", r ) );

        s.put("", parent_pop_seq_count );
        gen_seq_log.push_back( std::make_pair("", s ) );

        t.put("", nBlocks );
        gen_size_log.push_back( std::make_pair("", t ) );

        f.put("", fit_timer.elapsed().count());
        gen_fit_log.push_back( std::make_pair("", f ) );

        m.put("", repro_timer.elapsed().count());
        gen_repro_log.push_back( std::make_pair("", m ));

        c.put("", reset_timer.elapsed().count());
        gen_reset_log.push_back( std::make_pair("", c));
    }

    sim_time.stop();

    timer finalize_timer;
    unsigned int nUniqueInd = 0;
    for( unsigned int i = 1; i < parent->size(); ++i ) {
        bool no_dup = true;
        for( unsigned int j = 0; no_dup && j < i; ++j ) {
            no_dup = !((*parent)[i]->getProperties()->getGamete(0) == (*parent)[j]->getProperties()->getGamete(0) && (*parent)[i]->getProperties()->getGamete(1) == (*parent)[j]->getProperties()->getGamete(1))
                     && !((*parent)[i]->getProperties()->getGamete(0) == (*parent)[j]->getProperties()->getGamete(1) && (*parent)[i]->getProperties()->getGamete(1) == (*parent)[j]->getProperties()->getGamete(0));

        }
        if( no_dup ) {
            ++nUniqueInd;
        }
    }

    double nSymbols = 0;
    size_t nMaxSymbols = 0, nMinSymbols = -1;
    unsigned int n = 0;

    typedef std::vector< std::pair< locus_bitset::active_iterator, size_t > > count_vector_type;
    count_vector_type g_counts;

    std::vector< size_t > frequencies( alphabet_type::getInstance()->database_size(), 0);

    for( typename locus_bitset::active_iterator it = locus_bitset::active_begin(); it != locus_bitset::active_end(); it++ ) {
        (*it)->updateFrequencies( frequencies );
        size_t s = (*it)->size();
        nSymbols += (double)s;

        if( nMaxSymbols < s ) {
            nMaxSymbols = s;
        }

        if( nMinSymbols > s ) {
            nMinSymbols = s;
        }

        n += (*it)->copies();

        bool isUnique = true;
        for( count_vector_type::iterator it2 = g_counts.begin(); it2 != g_counts.end(); it2++ ) {
            if( (*it) == *(it2->first) ) {
                it2->second += (*it)->copies();
                isUnique = false;
                break;
            }
        }

        if( isUnique ) {
            g_counts.push_back(std::make_pair( it, (*it)->copies() ));
        }
    }

    double nGametes = (double)g_counts.size();

    boost::property_tree::ptree mut_stats;
    logMutationStats( mut_stats, frequencies, alphabet_type::getInstance() );

    n += (gamete_type::EMPTY.copies() - 1);
    alphabet_type::getInstance()->logState( log );

    log.put( "population.regions.count", n );
    log.put( "population.regions.unique", nGametes );
    log.put( "population.regions.active", locus_bitset::activeCount() );
    log.put( "population.regions.reference", gamete_type::EMPTY.copies() - 1);

    log.put( "population.individual.count", parent->size() );
    log.put( "population.individual.unique", nUniqueInd );
    log.put( "population.individual.selfing", nSelfing );

    log.put( "population.mutations.count", alphabet_type::getInstance()->active_count());
    log.put( "population.mutations.lost_or_fixed", alphabet_type::getInstance()->fixed_lost_count());
    log.put( "population.mutations.single_region.max", nMaxSymbols );
    log.put( "population.mutations.single_region.min", nMinSymbols );
    log.put( "population.mutations.single_region.average", nSymbols/nGametes );

    log.add_child( "population.mutations", mut_stats );

    log.put( "simulation.events.recombination", repro.getRecombinationEvents() );
    log.put( "simulation.events.mutation", repro.getMutationEvents() );

    log.put( "symbol_database.count", alphabet_type::getInstance()->database_size() );
    log.put( "symbol_database.symbol_size", sizeof( alphabet_type::symbol_type));
    log.put( "symbol_database.symbol_per_block", locus_bitset::bitset_type::bits_per_block );
    log.put( "symbol_database.max_block_per_region", (alphabet_type::getInstance()->database_size() / locus_bitset::bitset_type::bits_per_block) + 1 );
    log.put( "simulation.profile.memory_usage.total_blocks", nBlocks );
    log.put( "simulation.profile.memory_usage.blocks_per_region", nBlocks/nGametes );


#ifdef LOGGING
    boost::property_tree::write_json( std::cerr, global_log );
#endif

    while( !population.empty() ) {
        individual_type * ind = population.back();
        population.pop_back();

        delete ind;
    }

    while( !buffer.empty() ) {
        individual_type * ind = buffer.back();
        buffer.pop_back();
        delete ind;
    }

    finalize_timer.stop();
    runtime_timer.stop();

    log.put( "simulation.performance.initialize.nanosecond", init_time.elapsed().count() );
    log.put( "simulation.performance.simulate.nanoseconds", sim_time.elapsed().count() );
    log.put( "simulation.performance.finalize.nanoseconds", finalize_timer.elapsed().count() );
    log.put( "simulation.performance.runtime.nanoseconds", runtime_timer.elapsed().count() );
    log.add_child( "simulation.performance.data.generations", gen_log );
    log.add_child( "simulation.performance.data.runtimes", gen_run_log );
    log.add_child( "simulation.performance.data.table_size", gen_size_log );
    log.add_child( "simulation.performance.data.segregation_sites", gen_alpha_log );
    log.add_child( "simulation.performance.data.sequence_count", gen_seq_log );
    log.add_child( "simulation.performance.data.fitness", gen_fit_log);
    log.add_child( "simulation.performance.data.reproduction", gen_repro_log);
    log.add_child( "simulation.performance.data.reset", gen_reset_log);


    oss.str("");
    oss.clear();

    oss << vm[ OUTPUT_K ].as< string >() << "_bench_all_neutral.json";
    boost::property_tree::write_json( oss.str(), log );

    delete [] fitness;

    return 0;
}
