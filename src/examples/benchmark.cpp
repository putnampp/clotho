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

void write_log( boost::property_tree::ptree & log, const simulation_config & sim );

int main( int argc, char ** argv ) {
    simulation_config cfg;
    int res = parse_commandline(argc, argv, cfg);
    if( res ) {
        return res;
    }

    log_type    _log;
    rng_type    _rng;
    _rng.seed( cfg.seed );

    allele_set_type alleles;
    population_type pop, buffer;

    normal_dist_type ndist;
    allele_generator agen( ndist );

    predicate_type pred( &alleles );

    for( unsigned int i = 0; i < cfg.nPop; ++i ) {
        sequence_pointer r0 = alleles.empty_set();
        sequence_pointer r1 = alleles.empty_set();

        pop.push_back( std::make_pair( r0, r1 ) );
        buffer.push_back( std::make_pair(r0, r1) );
    }

    write_log( _log, cfg );   
    return res;
}

void write_log( boost::property_tree::ptree & log, const simulation_config & cfg ) {
    if( cfg.out_path.empty() ) {
        add_config( log, cfg );

        boost::property_tree::write_json(std::cout, log );
    } else {
        boost::property_tree::ptree _c;
        add_config( _c, cfg);

        string cpath = cfg.out_path + "_config.json";
        boost::property_tree::write_json( cpath, _c );

        cpath = cfg.out_path + "_log.json";
        boost::property_tree::write_json( cpath, log );
    }
}
