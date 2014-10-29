#include "config.h"

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <iterator>

#include "common_commandline.h"

#include "basic_allele.h"

#include "clotho/utility/timer.hpp"

#include "clotho/powerset/variable_subset.hpp"
#include "clotho/powerset/powerset_no_dup_pred.hpp"
#include "clotho/powerset/variable_subset_recombination.hpp"
#include "clotho/mutation/infinite_site_pred.hpp"
#include "clotho/mutation/element_generator.hpp"

#include "clotho/classifiers/region_classifier.hpp"

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

typedef clotho::classifiers::region_classifier< allele_type >               classifier_type;
typedef typename classifier_type::region_upper_bounds                       recombination_points;
typedef std::vector< allele_type >                                          allele_group;

typedef clotho::recombine::recombination< sequence_type, classifier_type >    recombination_type;

typedef std::pair< sequence_pointer, sequence_pointer >                     individual_type;
typedef std::vector< individual_type >                                      population_type;

typedef boost::random::poisson_distribution<unsigned int, double>           poisson_dist_type;
typedef boost::random::uniform_01< key_type >                               normal_dist_type;
typedef boost::random::uniform_int_distribution< unsigned int >             uniform_dist_type;

typedef clotho::mutations::element_generator< allele_type, normal_dist_type >   allele_generator;
typedef clotho::mutations::infinite_site_pred< allele_set_type >::type          predicate_type;

struct config_wrapper : public simulation_config {
    log_type    m_log;
    rng_type    m_rng;

    allele_set_type  alleles;
    allele_generator agen;
    predicate_type   pred;

    config_wrapper( const simulation_config & cfg ) :
        simulation_config( cfg )
        , m_log()
        , m_rng( seed )
        , alleles()
        , agen()
        , pred(&alleles)
    {}


    void reset_alleles() {
        alleles.clear();
    }

    virtual ~config_wrapper() { }
};

void write_log( config_wrapper & sim, unsigned int log_idx );
void initializePopulation(config_wrapper & cfg, population_type & pop, population_type & buffer, boost::property_tree::ptree & log );
void simulate( config_wrapper & cfg, population_type & pop, population_type & buffer, boost::property_tree::ptree & _log );
void resetPopulation( config_wrapper & cfg, population_type * p, boost::property_tree::ptree & _log );
void fitnessOfPopulation( config_wrapper & cfg, population_type * p, boost::property_tree::ptree & _log );
void reproducePopulation( config_wrapper & cfg, population_type * p, population_type * c, boost::property_tree::ptree & _log );

void statsPopulation( config_wrapper & cfg, population_type * p, boost::property_tree::ptree & _log );

sequence_pointer reproduce( config_wrapper & cfg, sequence_pointer s0, sequence_pointer s1, boost::property_tree::ptree & _log);
sequence_pointer recombine( config_wrapper & cfg, sequence_pointer base, sequence_pointer alt, unsigned int nEvents, bool copy = false);
void mutate( config_wrapper & cfg, sequence_pointer seq, unsigned int nEvents );


int main( int argc, char ** argv ) {
    simulation_config cmd;
    int res = parse_commandline(argc, argv, cmd);
    if( res ) {
        return res;
    }

    for( unsigned int i = 0; i < cmd.nRep; ++i ) {
        config_wrapper cfg(cmd);
        population_type pop, buffer;

        initializePopulation(cfg, pop, buffer, cfg.m_log );
        simulate( cfg, pop, buffer, cfg.m_log );

        write_log( cfg, i );
    }

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

    std::ostringstream oss;

    unsigned int log_period = ((cfg.nGen < cfg.log_period) ? cfg.nGen : cfg.log_period);

    for( unsigned int i = 0; i < cfg.nGen; ++i ) {
        boost::property_tree::ptree fitness_log;
        fitnessOfPopulation( cfg, parent, fitness_log );

        boost::property_tree::ptree repro_log;
        reproducePopulation(cfg, parent, child, repro_log );

        boost::property_tree::ptree reset_log;
        resetPopulation( cfg, parent, reset_log );
        std::swap( parent, child );

        cfg.alleles.pruneSpace();

        if( !(--log_period) ) {
            boost::property_tree::ptree stat_log;
            statsPopulation(cfg, parent, stat_log );

            if( !stat_log.empty() ) {
                oss.str("");
                oss.clear();
                oss << "data." << i;
                _log.add_child( oss.str(), stat_log );
            }

            log_period = ((i + cfg.log_period < cfg.nGen) ? cfg.log_period : (cfg.nGen - i - 1) );
        }
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
    poisson_dist_type mu_dist( cfg.mu ), rho_dist( cfg.rho );
    normal_dist_type  ndist;

    unsigned int nMut = mu_dist( cfg.m_rng );
    unsigned int nRec = rho_dist( cfg.m_rng );

    if( ndist( cfg.m_rng ) >= 0.5 ) {
        std::swap( s0, s1 );
    }
    
    sequence_pointer res;

    if( nRec ) {
        res = recombine( cfg, s0, s1, nRec, (nMut == 0) );
        if( nMut ) {
            assert( res );
            mutate( cfg, res, nMut );
        }
    } else if( nMut ) {
        // mutation only
        res = ((s0) ? s0->clone() : cfg.alleles.create_subset() );
        mutate(cfg, res, nMut );
    } else if( s0 ){
        // no recombination or mutation events
        s0->copy();
        res = s0;
    } 
//  else {
//      // implies that no recombination or mutations occurred and
//      // and that s0 was NULL
//      // therefore res should remain NULL as it was initialized
//    }

    return res;
}

sequence_pointer recombine( config_wrapper & cfg, sequence_pointer base, sequence_pointer alt, unsigned int nEvents, bool should_copy) {
    if( base == alt ) {
        if( !base ) {
            return cfg.alleles.create_subset();
        } else if( should_copy ) {
            base->copy();
            return base;
        } else {
            sequence_pointer res = base->clone();
            return res;
        }
    }

    recombination_points pts;
    allele_generator gen;

    gen( cfg.m_rng, nEvents, std::back_inserter(pts) );
    assert( pts.size() == nEvents );

    classifier_type cfier( pts );

    recombination_type rec;
    rec( base, alt, cfier );

    sequence_pointer res;
    if( !rec.isEmpty() ) {
        if( !should_copy ) {
            res = cfg.alleles.create_subset( *rec.getResultSequence() );
        } else if( rec.isMatchBase() ) {
            base->copy();
            res = base;
        } else if( rec.isMatchAlt() ) {
            alt->copy();
            res = base;
        } else {
            res = cfg.alleles.create_subset( *rec.getResultSequence() );
        }
    } else {
        res = cfg.alleles.create_subset();
    }

    return res;
}

void mutate( config_wrapper & cfg, sequence_pointer seq, unsigned int nEvents ) {
    allele_group pts;
    allele_generator gen;

    gen(cfg.m_rng, nEvents, &cfg.pred, std::back_inserter(pts) );

    assert( pts.size() == nEvents );

    while( !pts.empty() ) {
        seq->addElement( pts.back() );
        pts.pop_back();
    }
}

void statsPopulation( config_wrapper & cfg, population_type * p, boost::property_tree::ptree & _log ) {

    typedef std::map< unsigned int, unsigned int >  count_map;
    typedef count_map::iterator                     count_iterator;

    count_map allele_counts;

    typedef std::map< sequence_pointer, unsigned int >  sequence_map;
    typedef typename sequence_map::iterator             sequence_iterator;
    sequence_map seq_ref_counts;

    unsigned int nNull = 0, nSeq = 0;
    for( unsigned int i = 0; i < p->size(); ++i ) {
        if( !p->at(i).first ) {
            ++nNull;
        } else {
            assert( cfg.alleles.isFamilyMember( p->at(i).first ) );
            sequence_iterator sit = seq_ref_counts.find( p->at(i).first );

            if( sit != seq_ref_counts.end() ) {
                sit->second++;
            } else {
                seq_ref_counts.insert( std::make_pair( p->at(i).first, 1 ) );
            }
            ++nSeq;
        }

        if( !p->at(i).second ) {
            ++nNull;
        } else {
            assert( cfg.alleles.isFamilyMember( p->at(i).second ) );
            sequence_iterator sit = seq_ref_counts.find( p->at(i).second );

            if( sit != seq_ref_counts.end() ) {
                sit->second++;
            } else {
                seq_ref_counts.insert( std::make_pair( p->at(i).second, 1 ) );
            }
            ++nSeq;
        }
    }

    if(nNull) {
        allele_counts.insert( std::make_pair(0, nNull) );
    }

    unsigned int nExpSeq = nNull;
    double nAlleles = 0;

    bool usual_average = (nSeq % 2 == 0);
    unsigned int median_idx = nSeq / 2, med_allele = 0, prev_allele = 0;
    unsigned int mode_allele = 0, mode_seq = 0;
    for( sequence_iterator sit = seq_ref_counts.begin(); sit != seq_ref_counts.end(); sit++ ) {
        unsigned int acnt = sit->first->count();
        unsigned int e_rcnt = sit->second;

        count_iterator a_it = allele_counts.find( acnt );
        if( a_it == allele_counts.end() ) {
            allele_counts.insert( std::make_pair( acnt, e_rcnt ) );
        } else {
            a_it->second += e_rcnt;
        }
        nAlleles += (acnt * e_rcnt);

        unsigned int ubound = nExpSeq + e_rcnt;
        if(nExpSeq <= median_idx && median_idx < ubound) {
            if( !usual_average || ( nExpSeq <= median_idx - 1) ) {
                med_allele = acnt;
            } else {
                med_allele = (acnt + prev_allele) / 2;
            }
        } else {
            prev_allele = acnt;
        }

        if( e_rcnt > mode_seq ) {
            mode_allele = acnt;
            mode_seq = e_rcnt;
        }

        nExpSeq = ubound;
    }

    double ave_alleles_per_sequence = (nAlleles / (double) nExpSeq);

    // for validation purposes
    allele_set_type::cfamily_iterator f_it = cfg.alleles.family_begin(), f_end = cfg.alleles.family_end();
    unsigned int nNotInPop = 0;
    while( f_it != f_end ) {
        if( seq_ref_counts.find( *f_it++ ) == seq_ref_counts.end()) {
            ++nNotInPop;
        }
    }

    assert( nNotInPop == 0 );

    boost::property_tree::ptree all_dist;
    for( count_iterator c_it = allele_counts.begin(); c_it != allele_counts.end(); ++c_it) {
        boost::property_tree::ptree a,b,c;
        a.put( "", c_it->first );
        b.put( "", c_it->second );

        c.push_back (std::make_pair( "", a ) );
        c.push_back (std::make_pair( "", b ) );
        all_dist.push_back( std::make_pair("", c ) );
    }

    _log.put( "stats.population.family_size", cfg.alleles.family_size() );

    if( allele_counts.find(0) != allele_counts.end())
        _log.put( "stats.population.reference_sequences", allele_counts[0] );
    else
        _log.put( "stats.population.references_sequences", 0 );
    _log.put( "stats.population.total_sequences", nExpSeq );
    _log.put( "stats.sequences.comments", "Allele Count Distribution Format: [ allele_count, seq_ref_counts ]" );
    _log.add_child( "stats.sequences.allele_count_distribution", all_dist );
    _log.put( "stats.sequences.alleles_per.mean", ave_alleles_per_sequence);
    _log.put( "stats.sequences.alleles_per.min", allele_counts.begin()->first );
    _log.put( "stats.sequences.alleles_per.max", allele_counts.rbegin()->first );
    _log.put( "stats.sequences.alleles_per.median", med_allele );
    _log.put( "stats.sequences.alleles_per.mode", mode_allele );

    _log.put( "stats.alleles.variable_count", cfg.alleles.variable_size() );
    _log.put( "stats.alleles.fixed_count", cfg.alleles.fixed_size() );
}

void write_log( config_wrapper & cfg, unsigned int log_idx ) {
    if( cfg.out_path.empty() ) {
        add_config( cfg.m_log, cfg );

        boost::property_tree::write_json(std::cout, cfg.m_log );
    } else {
        boost::property_tree::ptree _c;
        add_config( _c, cfg);

        std::ostringstream oss;
        oss << cfg.out_path << "_" << log_idx;
        boost::property_tree::write_json( oss.str() + "_config.json", _c );
        boost::property_tree::write_json( oss.str() + "_log.json", cfg.m_log );
    }
}
