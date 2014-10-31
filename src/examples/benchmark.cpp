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
#include "clotho/powerset/variable_subset_fitness.hpp"
#include "clotho/mutation/infinite_site_pred.hpp"
#include "clotho/mutation/element_generator.hpp"

#include "clotho/classifiers/region_classifier.hpp"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

//typedef boost::property_tree::ptree                                         log_type;
typedef boost::random::mt19937                                              rng_type;
typedef clotho::utility::timer                                              timer_type;

typedef basic_allele                                                        allele_type;
typedef clotho::powersets::variable_subset< allele_type, unsigned long >    sequence_type;
typedef sequence_type::pointer                                              sequence_pointer;
typedef sequence_type::powerset_type                                        allele_set_type;

typedef typename allele_set_type::element_keyer_type                        element_keyer;
typedef typename element_keyer::key_type                                    key_type;

typedef clotho::classifiers::region_classifier< allele_type >               classifier_type;
typedef typename classifier_type::region_upper_bounds                       recombination_points;
typedef std::vector< allele_type >                                          allele_group;

typedef clotho::fitness::no_fit                                             het_fit_type;
typedef clotho::fitness::no_fit                                             alt_hom_fit_type;
typedef clotho::fitness::no_fit                                             ref_hom_fit_type;

typedef clotho::recombine::recombination< sequence_type, classifier_type >    recombination_type;

typedef clotho::fitness::fitness< sequence_type, het_fit_type, alt_hom_fit_type, ref_hom_fit_type, double > fitness_type;
typedef typename fitness_type::result_type                                  fitness_result_type;
typedef std::vector< fitness_result_type >                                  population_fitness_type;

typedef std::pair< sequence_pointer, sequence_pointer >                     individual_type;
typedef std::vector< individual_type >                                      population_type;
typedef typename population_type::iterator                                  population_iterator;

typedef std::vector< double >                                               population_fitness_type;

typedef boost::random::poisson_distribution<unsigned int, double>           poisson_dist_type;
typedef boost::random::uniform_01< key_type >                               normal_dist_type;
typedef boost::random::uniform_int_distribution< unsigned int >             uniform_dist_type;

typedef clotho::mutations::element_generator< allele_type, normal_dist_type >   allele_generator;
typedef clotho::mutations::infinite_site_pred< allele_set_type >::type          predicate_type;

struct config_wrapper : public simulation_config {
    boost::property_tree::ptree    m_log;
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

template < class Value >
void add_value_array( boost::property_tree::ptree & array, const Value & t ) {
    boost::property_tree::ptree node;

    node.put("", t );
    array.push_back( std::make_pair("", node));
}

void add_value_array( boost::property_tree::ptree & array, const clotho::utility::timer & t );

void add_node( boost::property_tree::ptree & array, const string & path, const boost::property_tree::ptree & n );

void initializePopulation(config_wrapper & cfg, population_type & pop, population_type & buffer, boost::property_tree::ptree & log );
void simulate( config_wrapper & cfg, population_type & pop, population_type & buffer, boost::property_tree::ptree & _log );
void resetPopulation( config_wrapper & cfg, population_type * p, boost::property_tree::ptree & _log );
double fitnessOfPopulation( config_wrapper & cfg, population_type * p, std::vector< double > & fit, boost::property_tree::ptree & _log );
void reproducePopulation( config_wrapper & cfg, population_type * p, population_type * c, boost::property_tree::ptree & _log );

void statsPopulation( config_wrapper & cfg, population_type * p, boost::property_tree::ptree & _log );

sequence_pointer reproduce( config_wrapper & cfg, sequence_pointer s0, sequence_pointer s1, boost::property_tree::ptree & _log);
sequence_pointer recombine( config_wrapper & cfg, sequence_pointer base, sequence_pointer alt, unsigned int nEvents, bool copy = false);
void mutate( config_wrapper & cfg, sequence_pointer seq, unsigned int nEvents );

typedef std::map< sequence_pointer, unsigned int >  ref_map_type;
typedef typename ref_map_type::iterator             ref_map_iterator;
void buildPopulationRefMap( population_type * p, ref_map_type & m);

struct scale_stats {
    boost::property_tree::ptree family_size, alleles_size, max_alleles;
};

void computeScaleStats( config_wrapper & cfg, population_type * p, scale_stats & stats );
void add_node( boost::property_tree::ptree &, const string & path, const scale_stats & s );

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
    sequence_pointer s = cfg.alleles.create_subset();
    for( unsigned int i = 0; i < cfg.nPop; ++i ) {
        pop.push_back( std::make_pair( s, s ) );
        buffer.push_back( std::make_pair(s, s ) );
    }

    _log.put( "status", "success" );
}

void simulate( config_wrapper & cfg, population_type & pop, population_type & buffer, boost::property_tree::ptree & _log ) {
    population_type *parent = &pop, *child = &buffer;

    std::ostringstream oss;

    unsigned int log_period = ((cfg.nGen < cfg.log_period) ? cfg.nGen : cfg.log_period);

    typedef boost::property_tree::ptree log_type;
    log_type perf_run, perf_fit, perf_repro, perf_reset, perf_stats;
    log_type l_fit;

    scale_stats mem_stats;

    for( unsigned int i = 0; i < cfg.nGen; ++i ) {
        oss.str("");
        oss.clear();
        oss << "simulation." << i;

        timer_type run_t;
        timer_type fit_t;

        population_fitness_type fit;
        boost::property_tree::ptree fit_log;
        double pop_fit = fitnessOfPopulation( cfg, parent, fit, fit_log );
        double e_fit = ((double) parent->size() / pop_fit);
        fit_t.stop();

        timer_type repro_t;
        boost::property_tree::ptree repro_log;
        reproducePopulation(cfg, parent, child, repro_log );
        repro_t.stop();

        timer_type reset_t;
        boost::property_tree::ptree reset_log;
        resetPopulation( cfg, parent, reset_log );
        std::swap( parent, child );

        cfg.alleles.pruneSpace();

        reset_t.stop();

        timer_type stat_t;
        boost::property_tree::ptree stat_log;
        if( !(--log_period) ) {
            std::cerr << "computing population stats: " << i << std::endl;
            statsPopulation(cfg, parent, stat_log );
            stat_t.stop();
            log_period = ((i + cfg.log_period < cfg.nGen) ? cfg.log_period : (cfg.nGen - i - 1) );
        }
        stat_t.stop();
        run_t.stop();

        computeScaleStats( cfg, parent, mem_stats );

        add_value_array( l_fit, e_fit);

        add_value_array( perf_run, run_t );
        add_value_array( perf_fit, fit_t);
        add_value_array( perf_repro, repro_t);
        add_value_array( perf_reset, reset_t);
        add_value_array( perf_stats, stat_t);

        add_node( _log, oss.str() + ".fitness", fit_log);
        add_node( _log, oss.str() + ".reproduction", repro_log);
        add_node( _log, oss.str() + ".reset", reset_log);
        add_node( _log, oss.str() + ".stats", stat_log );
    }

    const double h = timer_type::hertz;
    const string sim_path = "simulation.performance.data";
    _log.put( sim_path + ".runtimes.scale", h );
    _log.add_child( sim_path + ".runtimes.total", perf_run );
    _log.add_child( sim_path + ".runtimes.fitness", perf_fit );
    _log.add_child( sim_path + ".runtimes.reproduction", perf_repro );
    _log.add_child( sim_path + ".runtimes.reset", perf_reset );
    _log.add_child( sim_path + ".runtimes.stats", perf_stats );

    _log.add_child("simulation.data.fitness", l_fit );

    add_node( _log, sim_path, mem_stats );
}

inline void add_value_array( boost::property_tree::ptree & array, const clotho::utility::timer & t ) {
    boost::property_tree::ptree node;

    node.put("", t.elapsed().count() );
    array.push_back( std::make_pair("", node));
}

inline void add_node( boost::property_tree::ptree & root, const string & path, const boost::property_tree::ptree & node ) {
    if( !node.empty() ) {
        root.add_child( path, node );
    }
}

void computeScaleStats( config_wrapper & cfg, population_type * p, scale_stats & stats ) {
    ref_map_type seq_ref_counts;

    size_t  fam_size = 0;
    
    allele_set_type::cfamily_iterator f_it = cfg.alleles.family_begin(), f_end = cfg.alleles.family_end();
    while( f_it != f_end ) {
        if( (*f_it) ) {
            fam_size += (*f_it)->num_blocks();
        }
        ++f_it;
    }

    size_t  max_allele_size = cfg.alleles.variable_allocated_size();
    size_t  allele_size = cfg.alleles.variable_size();

    add_value_array( stats.family_size, fam_size );
    add_value_array( stats.alleles_size, allele_size );
    add_value_array( stats.max_alleles, max_allele_size);
}

void add_node( boost::property_tree::ptree & r, const string & path, const scale_stats & s ) {
    add_node( r, path + ".memory.family_size", s.family_size );
    add_node( r, path + ".memory.alleles_size", s.alleles_size );
    add_node( r, path + ".memory.max_alleles", s.max_alleles );
}

void resetPopulation( config_wrapper & cfg, population_type * p, boost::property_tree::ptree & _log ) {
    population_type::iterator it = p->begin();
    while( it != p->end() ) {
        it->first.reset();
        it->second.reset();
        ++it;
    }
}

double fitnessOfPopulation( config_wrapper & cfg, population_type * p, std::vector< double > & fit, boost::property_tree::ptree & _log ) {
    population_iterator first = p->begin();

    double pop_fit = 0.;
    fitness_type    f;
    std::back_insert_iterator< std::vector< double > > output = std::back_inserter(fit);
    while( first != p->end() ) {
        double val = f(0., first->first, first->second);
        pop_fit += (val);
        (*output++) = val;
        ++first;
    }
    return pop_fit;
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

        if( !p0_log.empty() || !p1_log.empty() ) {
            std::ostringstream oss("child.");
            oss << (it - c->begin());
            add_node( _log, oss.str() + ".p0", p0_log);
            add_node( _log, oss.str() + ".p1", p1_log);
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
            res = base;
        } else if( rec.isMatchAlt() ) {
            res = alt;
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

void buildPopulationRefMap( population_type * p, ref_map_type & m ) {
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
}

void statsPopulation( config_wrapper & cfg, population_type * p, boost::property_tree::ptree & _log ) {

    typedef std::map< unsigned int, unsigned int >  count_map;
    typedef count_map::iterator                     count_iterator;

    count_map allele_counts;

    ref_map_type seq_ref_counts;
    buildPopulationRefMap( p, seq_ref_counts );

    unsigned int nExpSeq = 0;
    unsigned int nSeq = seq_ref_counts.size(); // - ((nExpSeq > 0)? 1 : 0); // remove null sequence
    double nAlleles = 0;

    bool usual_average = ( nSeq % 2 == 0);
    unsigned int median_idx = nSeq / 2, med_allele = 0, prev_allele = 0;
    unsigned int mode_allele = 0, mode_seq = 0;
    for( ref_map_iterator sit = seq_ref_counts.begin(); sit != seq_ref_counts.end(); sit++ ) {
        unsigned int acnt = (( sit->first ) ? sit->first->count() : 0 );
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

    _log.put( "stats.population.references_sequences", ((allele_counts.find(0) != allele_counts.end())? allele_counts[0] : 0 ) );
    
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
