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

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <iterator>

#include "common_commandline.h"

#include "basic_allele.h"

#include "clotho/utility/timer.hpp"

#include "clotho/utility/bit_block_iterator.hpp"
#include "clotho/powerset/variable_subset.hpp"
#include "clotho/powerset/powerset_no_dup_pred.hpp"
#include "clotho/powerset/variable_subset_recombination.hpp"
#include "clotho/powerset/variable_subset_fitness.hpp"
#include "clotho/mutation/infinite_site_pred.hpp"
#include "clotho/mutation/element_generator.hpp"

#include "clotho/classifiers/region_classifier.hpp"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/discrete_distribution.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

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

typedef clotho::fitness::fitness_method< double, clotho::fitness::multiplicative_heterozygous_tag > het_fit_type;
typedef clotho::fitness::fitness_method< double, clotho::fitness::multiplicative_homozygous_tag >    alt_hom_fit_type;
typedef clotho::fitness::no_fit                                             ref_hom_fit_type;

typedef clotho::recombine::recombination< sequence_type, classifier_type, clotho::recombine::inspection::tag::copy_matching_classify_mismatch, clotho::recombine::walker::tag::inline_dynamic_classify >    recombination_type;

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
typedef boost::random::discrete_distribution< unsigned int, double >        discrete_dist_type;

typedef clotho::mutations::infinite_site_pred< allele_set_type >::type          predicate_type;

struct config_wrapper : public simulation_config {
    boost::property_tree::ptree    m_log;
    rng_type    m_rng;

    allele_set_type  alleles;
    predicate_type   pred;

    config_wrapper( const simulation_config & cfg ) :
        simulation_config( cfg )
        , m_log()
        , m_rng( seed )
        , alleles()
        , pred(&alleles) {
    }


    void reset_alleles() {
        alleles.clear();
    }

    virtual ~config_wrapper() { }
};

template < class KeyDist >
struct custom_allele_gen {
    typedef KeyDist key_dist_type;

    key_dist_type m_kd;

    custom_allele_gen() {}

    custom_allele_gen( const key_dist_type & kd ) :
        m_kd(kd) {
    }

    custom_allele_gen( const custom_allele_gen< KeyDist > & cag ) :
        m_kd( cag.m_kd ) {
    }

    template < class URNG >
    basic_allele operator()( URNG & rng, double sel = 0.0, double dom = 1.0, bool neut = true ) {
#ifdef ALL_SELECTED_ALLELES
        neut = false;
#endif  // ALL_SELECTED_ALLELES
        return basic_allele( m_kd(rng), sel, dom, neut );
    }
};

//typedef clotho::mutations::element_generator< allele_type, normal_dist_type >   allele_generator;
typedef clotho::mutations::element_generator< allele_type, custom_allele_gen< normal_dist_type > >   allele_generator;

void write_log( config_wrapper & sim, unsigned int log_idx );

template < class Value >
void add_value_array( boost::property_tree::ptree & array, const Value & t ) {
    boost::property_tree::ptree node;

    node.put("", t );
    array.push_back( std::make_pair("", node));
}

template < class A, class B >
void add_value_array( boost::property_tree::ptree & array, const std::pair< A, B > & t ) {
    boost::property_tree::ptree a,b,c;
    a.put( "", t.first );
    b.put( "", t.second );

    c.push_back (std::make_pair( "", a ) );
    c.push_back (std::make_pair( "", b ) );
    array.push_back( std::make_pair("", c ) );
}

void add_value_array( boost::property_tree::ptree & array, const clotho::utility::timer & t );

void add_node( boost::property_tree::ptree & array, const string & path, const boost::property_tree::ptree & n );
void add_node( boost::property_tree::ptree & array, const string & path, const individual_type & ind  );

void initializePopulation(config_wrapper & cfg, population_type & pop, population_type & buffer, boost::property_tree::ptree & log );
void simulate( config_wrapper & cfg, population_type & pop, population_type & buffer, boost::property_tree::ptree & _log );
void resetPopulation( config_wrapper & cfg, population_type * p, boost::property_tree::ptree & _log );
double fitnessOfPopulation( config_wrapper & cfg, population_type * p, std::vector< double > & fit, boost::property_tree::ptree & _log );

void statsPopulation( config_wrapper & cfg, population_type * p, boost::property_tree::ptree & _log );

sequence_pointer reproduce( config_wrapper & cfg, sequence_pointer s0, sequence_pointer s1, boost::property_tree::ptree & _log);
sequence_pointer recombine( config_wrapper & cfg, sequence_pointer base, sequence_pointer alt, unsigned int nEvents, bool copy = false);
void mutate( config_wrapper & cfg, sequence_pointer seq, unsigned int nEvents );

typedef std::map< sequence_pointer, unsigned int >  ref_map_type;
typedef typename ref_map_type::iterator             ref_map_iterator;
typedef std::vector< unsigned int >                 allele_dist_type;
void buildPopulationRefMap( population_type * p, ref_map_type & m, allele_dist_type & a);

struct scale_stats {
    boost::property_tree::ptree family_size, alleles_size, max_alleles;
};

void computeScaleStats( config_wrapper & cfg, population_type * p, scale_stats & stats );
void add_node( boost::property_tree::ptree &, const string & path, const scale_stats & s );

template < class PopDist >
void reproducePopulation( config_wrapper & cfg, population_type * p, population_type * c, const PopDist & dist, boost::property_tree::ptree & _log ) {
    population_type::iterator it = c->begin();
    while( it != c->end() ) {
        unsigned int idx = dist( cfg.m_rng ), idx2 = dist( cfg.m_rng );

        population_type::iterator p0 = p->begin() + idx, p1 = p->begin() + idx2;

        boost::property_tree::ptree p0_log;
        it->first = reproduce( cfg, p0->first, p0->second, p0_log );

        boost::property_tree::ptree p1_log;
        it->second = reproduce( cfg, p1->first, p1->second, p1_log );

        if( !p0_log.empty() || !p1_log.empty() ) {
            std::ostringstream oss;
            oss << "child." << (it - c->begin());
            add_node( _log, oss.str() + ".p0", p0_log);
            add_node( _log, oss.str() + ".p1", p1_log);
        }

        ++it;
    }
}

int main( int argc, char ** argv ) {
    simulation_config cmd;
    int res = parse_commandline(argc, argv, cmd);
    if( res ) {
        return res;
    }

    for( unsigned int i = 0; i < cmd.nRep; ++i ) {
        config_wrapper cfg(cmd);

        if( cfg.nGen > 0 ) {
            population_type pop, buffer;

            initializePopulation(cfg, pop, buffer, cfg.m_log );
            simulate( cfg, pop, buffer, cfg.m_log );
        }

        write_log( cfg, i );
    }

    return res;
}

void initializePopulation( config_wrapper & cfg, population_type & pop, population_type & buffer, boost::property_tree::ptree & _log ) {
    for( unsigned int i = 0; i < cfg.nPop; ++i ) {
        pop.push_back( std::make_pair( cfg.alleles.create_subset() , cfg.alleles.create_subset()  ) );
        buffer.push_back( std::make_pair( cfg.alleles.create_subset(), cfg.alleles.create_subset() ) );
    }

    _log.put( "status", "success" );
}

void simulate( config_wrapper & cfg, population_type & pop, population_type & buffer, boost::property_tree::ptree & _log ) {
    population_type *parent = &pop, *child = &buffer;

    std::ostringstream oss;

    unsigned int log_period = ((cfg.nGen < cfg.log_period) ? cfg.nGen : cfg.log_period);

    typedef boost::property_tree::ptree log_type;
    log_type perf_run, perf_fit, perf_repro, perf_reset, perf_stats;
    log_type l_fit, l_gen, l_fam;

    scale_stats mem_stats;

    for( unsigned int i = 0; i < cfg.nGen; ++i ) {
        oss.str("");
        oss.clear();
        oss << "simulation." << i;

        timer_type run_t;
        timer_type fit_t;

        population_fitness_type fit;
        boost::property_tree::ptree fit_log;
        double e_fit = 1.;
        double pop_fit = fitnessOfPopulation( cfg, parent, fit, fit_log );
        e_fit = ((double) parent->size() / pop_fit);
        fit_t.stop();

        timer_type repro_t;
        boost::property_tree::ptree repro_log;

        discrete_dist_type dist( fit.begin(), fit.end() );

        assert( dist.probabilities().size() == fit.size());
        reproducePopulation(cfg, parent, child, dist, repro_log );
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
        add_value_array( l_gen, i);
        add_value_array( l_fam, cfg.alleles.family_size() );

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
    _log.put( sim_path + ".runtime.scale", h );
    _log.add_child( sim_path + ".generations", l_gen );
    _log.add_child( sim_path + ".runtime.total", perf_run );
    _log.add_child( sim_path + ".runtime.fitness", perf_fit );
    _log.add_child( sim_path + ".runtime.reproduction", perf_repro );
    _log.add_child( sim_path + ".runtime.reset", perf_reset );
    _log.add_child( sim_path + ".runtime.stats", perf_stats );
    _log.add_child( sim_path + ".family_size", l_fam );

    _log.add_child(sim_path + ".e_fitness", l_fit );

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
    add_node( r, path + ".memory.table_size", s.family_size );
    add_node( r, path + ".memory.alleles_size", s.alleles_size );
    add_node( r, path + ".memory.max_alleles", s.max_alleles );
}

void add_node( boost::property_tree::ptree & r, const string & path, const individual_type & ind ) {
    std::ostringstream oss;

    oss << ind.first << ": " << (*ind.first);
    r.put( path + ".first", oss.str());

    oss.str("");
    oss.clear();

    oss << ind.second << ": " << (*ind.second);
    r.put( path + ".second", oss.str() );
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

    // pre-compute which alleles are selected within the population
    //
    // identify all variable (non-free) indicies (elements)
    typename allele_set_type::bitset_type selected( cfg.alleles.free_begin(), cfg.alleles.free_end());
    selected.flip();

    // foreach variable element
    size_t idx = selected.find_first();
    while( idx != allele_set_type::bitset_type::npos ) {
        typename allele_set_type::cvariable_iterator vit =  cfg.alleles.variable_begin() + idx;

        // determine if neutral
        if( vit->isNeutral() ) {
            // it is neutral therefore should be ignored in fitness calculation
            selected.reset(idx);
        }
        idx = selected.find_next(idx);
    }

    size_t nSelected = selected.count();

    std::back_insert_iterator< std::vector< double > > output = std::back_inserter(fit);

    if( nSelected == 0 ) {
        population_iterator first = p->begin();
        while( first != p->end() ) {
            (*output++) = 1.;
            ++first;
        }
        return p->size();
    }

    population_iterator first = p->begin();

    double pop_fit = 0.;
    fitness_type    f;

    ref_map_type m_refs;
    while( first != p->end() ) {
        double val = f(1., first->first, first->second, selected.m_bits.begin(), selected.m_bits.end());
        pop_fit += (val);
        (*output++) = val;
        ++first;
    }
    return pop_fit;
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
            assert( res && s0 != res && s1 != res );
            mutate( cfg, res, nMut );
        }
    } else if( nMut ) {
        // mutation only
        res = ((s0) ? s0->clone() : cfg.alleles.create_subset() );
        assert( res != s0 );
        mutate(cfg, res, nMut );
    } else if( s0 ) {
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

    // without predicate, therefore only random points are created
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
        typename sequence_type::index_type idx = rit->first->find_first();
        while( idx.second !=  sequence_type::npos ) {
            a[ idx.first ] += n;
            idx = rit->first->find_next( idx );
        }
    }
}

void statsPopulation( config_wrapper & cfg, population_type * p, boost::property_tree::ptree & _log ) {

    typedef std::map< unsigned int, unsigned int >  count_map;
    typedef count_map::iterator                     count_iterator;

    count_map allele_counts;

    ref_map_type seq_ref_counts;
    allele_dist_type allele_dist( cfg.alleles.variable_allocated_size(), 0);
    buildPopulationRefMap( p, seq_ref_counts, allele_dist );

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

    boost::property_tree::ptree all_dist, all_freq;
    for( count_iterator c_it = allele_counts.begin(); c_it != allele_counts.end(); ++c_it) {
        add_value_array( all_dist, *c_it );
    }

    typename allele_set_type::cvariable_iterator e_it = cfg.alleles.variable_begin();
    for( allele_dist_type::iterator a_it = allele_dist.begin(); a_it != allele_dist.end(); ++a_it, ++e_it ) {
        if( (*a_it) > 0 ) {
            std::ostringstream oss;
            oss << *e_it;
            add_value_array( all_freq, std::make_pair(oss.str(), (*a_it)));
        }
    }

    _log.put( "population.family_size", cfg.alleles.family_size() );

    _log.put( "population.reference_sequences", ((allele_counts.find(0) != allele_counts.end())? allele_counts[0] : 0 ) );

    _log.put( "population.total_sequences", nExpSeq );
    _log.put( "sequences.allele_count_distribution_format", "[alleles per sequence, number of sequences in population]" );
    _log.add_child( "sequences.allele_count_distribution", all_dist );
    _log.put( "sequences.alleles_per.mean", ave_alleles_per_sequence);
    _log.put( "sequences.alleles_per.min", allele_counts.begin()->first );
    _log.put( "sequences.alleles_per.max", allele_counts.rbegin()->first );
    _log.put( "sequences.alleles_per.median", med_allele );
    _log.put( "sequences.alleles_per.mode", mode_allele );

    _log.put( "alleles.variable_count", cfg.alleles.variable_size() );
    _log.put( "alleles.fixed_count", cfg.alleles.fixed_size() );
    _log.put( "alleles.free_size", cfg.alleles.free_size() );

    _log.put( "alleles.frequencies_format", "[allele, count in population]" );
    _log.add_child( "alleles.frequencies", all_freq );
}

void write_log( config_wrapper & cfg, unsigned int log_idx ) {
    boost::property_tree::ptree clog;
    add_config( clog, cfg );
    cfg.m_log.add_child( CONFIG_BLOCK_K, clog );

    if( cfg.out_path.empty() ) {
        boost::property_tree::write_json(std::cout, cfg.m_log );
    } else {
        std::ostringstream oss;
        oss << cfg.out_path << "_" << log_idx;
        boost::property_tree::write_json( oss.str() + "_log.json", cfg.m_log );
    }
}
