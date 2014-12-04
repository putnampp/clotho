#ifndef SIMULATE_ENGINE_HPP_
#define SIMULATE_ENGINE_HPP_

#include <algorithm>
#include <boost/property_tree/ptree.hpp>

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

#include "infinite_site.hpp"

#include "clotho/genetics/recombiner.hpp"

#include "clotho/powerset/powerset_no_dup_pred.hpp"
#include "clotho/powerset/variable_subset_recombination.hpp"
#include "clotho/powerset/variable_subset_fitness.hpp"

#include "clotho/utility/parameter_space.hpp"

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
    typedef typename base_type::population_type population_type;

    typedef clotho::utility::random_generator< rng_type, infinite_site< sequence_type > >       mutation_generator_type;
    typedef sequence_generator< sequence_pointer >                                              sequence_generator_type;
    typedef sequence_mutator< sequence_type, mutation_generator_type >                          sequence_mutator_type;
    typedef clotho::utility::random_generator< rng_type, sequence_mutator_type >                sequence_mutator_generator;
    typedef individual_initializer< individual_type, sequence_generator< sequence_pointer > >   individual_initializer_type;
    typedef individual_selector< rng_type >             individual_selector_type;
    typedef individual_reproduction< individual_type
                , sequence_mutator_generator
                , recombination_method_generator > individual_reproduction_type;
    typedef individual_generator< population_type, individual_selector_type, individual_reproduction_type >     individual_generator_type;

    typedef individual_resetter< individual_type >  individual_resetter_type;

    // fitness typedefs
    typedef clotho::fitness::fitness_method< double, clotho::fitness::multiplicative_heterozygous_tag > het_fit_type;
    typedef clotho::fitness::fitness_method< double, clotho::fitness::multiplicative_homozygous_tag >    alt_hom_fit_type;
    typedef clotho::fitness::no_fit                                             ref_hom_fit_type;
    typedef clotho::fitness::fitness< sequence_type, het_fit_type, alt_hom_fit_type, ref_hom_fit_type, double > fitness_type;
    typedef typename fitness_type::result_type                                  fitness_result_type;
    typedef std::vector< fitness_result_type >                                  population_fitness_type;
    typedef individual_fitness< fitness_type >                                  fitness_operator;

    simulate_engine( boost::property_tree::ptree & config ) :
        m_rng()
        , m_founder_size( DEFAULT_POPULATION_SIZE )
        , m_seq_mut_gen( m_rng, config )
        , m_rec_met_gen( m_rng, config )
        , m_repro( m_seq_mut_gen, m_rec_met_gen )
        , m_parent( &m_pop )
        , m_child( &m_buffer )
    {
        parseConfig( config );
        initialize();
    }

    void simulate() {
        simulate( m_founder_size );
    }

    void simulate( unsigned int p_size ) {
        population_fitness_type pfit;
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

        individual_selector_type sel( m_rng, pfit.begin(), pfit.end() );
        individual_generator_type ind_gen( m_parent, sel, m_repro );

        std::generate_n( std::back_inserter( *m_child ), p_size, ind_gen );
    }

    void reset_parent() {
        m_parent->clear();

        std::swap( m_parent, m_child );

        m_alleles.pruneSpace(); 
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

protected:
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
    }

    void initialize( ) {
        m_pop.clear();
        m_buffer.clear();

        m_pop.reserve( m_founder_size );
        m_buffer.reserve( m_founder_size );

        sequence_generator_type sgen( m_alleles );
        individual_initializer_type igen( sgen );
        std::generate_n( std::back_inserter( m_pop ), m_founder_size, igen );
    }

    rng_type        m_rng;
    unsigned int    m_founder_size;
    sequence_mutator_generator m_seq_mut_gen;
    recombination_method_generator m_rec_met_gen;

    individual_reproduction_type m_repro;

    log_type            m_log;
    population_type     m_pop, m_buffer;
    population_type     * m_parent, * m_child;

    allele_set_type m_alleles;
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
