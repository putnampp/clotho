#include "config.h"

#include <iostream>
#include <vector>

#include "common_commandline.h"

#include "basic_allele.h"

#include "clotho/utility/timer.hpp"

#include "clotho/powerset/variable_subset.hpp"
#include "clotho/powerset/powerset_no_dup_pred.hpp"
#include "clotho/mutation/infinite_site_pred.hpp"
#include "clotho/mutation/element_generator.hpp"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

typedef boost::property_tree::ptree                                         log_type;
typedef boost::random::mt19937                                              rng_type;

typedef basic_allele                                                        allele_type;
typedef clotho::powersets::variable_subset< allele_type, unsigned long >    sequence_type;
typedef sequence_type::pointer                                              sequence_pointer;
typedef sequence_type::powerset_type                                        allele_set_type;

typedef typename allele_set_type::element_keyer_type                        element_keyer;
typedef typename element_keyer::key_type                                    key_type;

typedef std::pair< sequence_pointer, sequence_pointer >                     individual_type;
typedef std::vector< individual_type >                                      population_type;

typedef boost::random::poisson_distribution<unsigned int, double>           poisson_dist_type;
typedef boost::random::uniform_01< key_type >                               normal_dist_type;
typedef boost::random::uniform_int_distribution< unsigned int >             uniform_dist_type;

typedef clotho::mutations::element_generator< basic_allele, normal_dist_type >  allele_generator;
typedef clotho::mutations::infinite_site_pred< allele_set_type >::type          predicate_type;

struct config_wrapper : public simulation_config {
    log_type    m_log;
    rng_type    m_rng;

    allele_generator agen;

    config_wrapper( const simulation_config & cfg ) :
        simulation_config( cfg )
        , m_log()
        , m_rng( seed )
        , agen()
    {}
};

void write_log( config_wrapper & sim );
void initializePopulation(config_wrapper & cfg, population_type & pop, population_type & buffer, boost::property_tree::ptree & log );
void simulate( config_wrapper & cfg, population_type & pop, population_type & buffer, boost::property_tree::ptree & _log );
void resetPopulation( config_wrapper & cfg, population_type * p, boost::property_tree::ptree & _log );
void fitnessOfPopulation( config_wrapper & cfg, population_type * p, boost::property_tree::ptree & _log );
void reproducePopulation( config_wrapper & cfg, population_type * p, population_type * c, boost::property_tree::ptree & _log );

sequence_pointer reproduce( config_wrapper & cfg, sequence_pointer s0, sequence_pointer s1, boost::property_tree::ptree & _log);


int main( int argc, char ** argv ) {
    simulation_config cmd;
    int res = parse_commandline(argc, argv, cmd);
    if( res ) {
        return res;
    }

    config_wrapper cfg(cmd);

    for( unsigned int i = 0; i < cfg.nRep; ++i ) {
        population_type pop, buffer;

        log_type init_log;
        initializePopulation(cfg, pop, buffer, init_log );
        
        log_type sim_log;
        simulate( cfg, pop, buffer, sim_log );

        std::ostringstream oss;
        oss << "repeat." << i;
        if( !init_log.empty() ) {
            cfg.m_log.add_child( oss.str() + ".init", init_log );
        }

        if( !sim_log.empty() ) {
            cfg.m_log.add_child( oss.str() + ".simulate", sim_log );
        }
    }

    write_log( cfg );   
    return res;
}

void initializePopulation( config_wrapper & cfg, population_type & pop, population_type & buffer, boost::property_tree::ptree & _log ) {
    for( unsigned int i = 0; i < cfg.nPop; ++i ) {
        pop.push_back( std::make_pair( sequence_pointer(), sequence_pointer() ) );
        buffer.push_back( std::make_pair(sequence_pointer(), sequence_pointer() ) );
    }

    _log.put( "status", "success" );
}

void simulate( config_wrapper & cfg, population_type & pop, population_type & buffer, boost::property_tree::ptree & _log ) {
    population_type *parent = &pop, *child = &buffer;

    allele_set_type alleles;
    predicate_type pred( &alleles );

    for( unsigned int i = 0; i < cfg.nGen; ++i ) {
        boost::property_tree::ptree fitness_log;
        fitnessOfPopulation( cfg, parent, fitness_log );

        boost::property_tree::ptree repro_log;
        reproducePopulation(cfg, parent, child, repro_log );

        boost::property_tree::ptree reset_log;
        resetPopulation( cfg, parent, reset_log );
        std::swap( parent, child );

        alleles.pruneSpace();
    }
}

void resetPopulation( config_wrapper & cfg, population_type * p, boost::property_tree::ptree & _log ) {
    population_type::iterator it = p->begin();
    while( it != p->end() ) {
        if( it->first )     it->first->release();
        it->first.reset();

        if( it->second )    it->second->release();
        it->second.reset();
        ++it;
    }
}

void fitnessOfPopulation( config_wrapper & cfg, population_type * p, boost::property_tree::ptree & _log ) {
}

void reproducePopulation( config_wrapper & cfg, population_type * p, population_type * c, boost::property_tree::ptree & _log ) {
    uniform_dist_type udist(0, p->size());

    population_type::iterator it = c->begin();
    while( it != c->end() ) {
        unsigned int idx = udist( cfg.m_rng ), idx2 = udist( cfg.m_rng );
        
        population_type::iterator p0 = p->begin() + idx, p1 = p->begin() + idx2;

        boost::property_tree::ptree p0_log;
        it->first = reproduce( cfg, p0->first, p0->second, p0_log );
        
        boost::property_tree::ptree p1_log;
        it->second = reproduce( cfg, p1->first, p1->second, p1_log );

        if( !p0_log.empty() ) {

        }

        if( !p1_log.empty() ) {

        }

        ++it;
    }
}

sequence_pointer reproduce( config_wrapper & cfg, sequence_pointer s0, sequence_pointer s1, boost::property_tree::ptree & _log) {
    if( s0 == s1 ) {
        if( s0 ) s0->copy();
        return s0;
    }
    poisson_dist_type mu_dist( cfg.mu ), rho_dist( cfg.rho );
    normal_dist_type  ndist;

    unsigned int nMut = mu_dist( cfg.m_rng );
    unsigned int nRec = rho_dist( cfg.m_rng );

    if( ndist( cfg.m_rng ) >= 0.5 ) {
        std::swap( s0, s1 );
    }
    
    if( nMut == 0 && nRec == 0 ) {
        // no recombination or mutation events
        s0->copy();
        return s0;
    }

    sequence_pointer res;

    if( nRec ) {
        if( nMut ) {
            res = recombine_and_mutate( s0, s1, nRec, nMut );
        } else {
            res = recombine( s0, s1, nRec );
        }
    } else {
        // mutation only
        res = s0->clone();        
    }

    return res;
}

void write_log( config_wrapper & cfg ) {
    if( cfg.out_path.empty() ) {
        add_config( cfg.m_log, cfg );

        boost::property_tree::write_json(std::cout, cfg.m_log );
    } else {
        boost::property_tree::ptree _c;
        add_config( _c, cfg);

        string cpath = cfg.out_path + "_config.json";
        boost::property_tree::write_json( cpath, _c );

        cpath = cfg.out_path + "_log.json";
        boost::property_tree::write_json( cpath, cfg.m_log );
    }
}
