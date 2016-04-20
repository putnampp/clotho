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
#include <boost/test/unit_test.hpp>

#include "../unittest_config.h"

#include "clotho/data_spaces/allele_space/allele_space.hpp"
#include "clotho/data_spaces/population_space/genetic_space.hpp"
#include "clotho/data_spaces/crossover/crossover.hpp"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/random/mersenne_twister.hpp>

using namespace clotho::genetics;

typedef boost::random::mt19937                              random_engine_type;
typedef double                                              position_type;
typedef double                                              weight_type;
typedef unsigned long long                                  block_type;
typedef qtl_allele_vectorized< position_type, weight_type > allele_type;
typedef genetic_space< allele_type, block_type >             genetic_space_type;

BOOST_AUTO_TEST_SUITE( test_data_space )

BOOST_AUTO_TEST_CASE( crossover_test ) {
    typedef Crossover< random_engine_type, genetic_space_type > crossover_type;

    typedef typename crossover_type::mate_pair_type mate_pair_type;

    boost::property_tree::ptree config;
    config.put( "recombination.rho", 1.0 );

    random_engine_type          rng(1234);

    genetic_space_type  parents, offspring;
    crossover_type      xover( &rng, config );

    boost::property_tree::write_json( std::cout, config );

    size_t exp_parents = 2, exp_alleles = 150, exp_offspring = 1;

    parents.grow( exp_parents, exp_alleles);

    parents.getSequenceSpace().clear();

    offspring.grow( exp_offspring, exp_alleles );
    offspring.getSequenceSpace().clear();

    mate_pair_type      mates;
    mates.push_back( std::make_pair( (size_t) 0, (size_t) 1) );

    // based upon random number generator initialization
    // and recombination rate configuration
    // an event will be generated when recombining parent genome 1
    // and the event will occur at 0.622109
    // any position greater than this should come from  sequence 3
    parents.getAlleleSpace().setPositionAt( 88, 0.75 );
    parents.getAlleleSpace().setPositionAt( 99, 0.85 );

    // embed some alleles in the parental genomes
    parents.getSequenceSpace().flip( 0, 28 );
    parents.getSequenceSpace().flip( 0, 110 );
    
    parents.getSequenceSpace().flip( 1, 28 );   // parent 0 homozygous at allele 28

    parents.getSequenceSpace().flip( 2, 77 );
    parents.getSequenceSpace().flip( 2, 99 );   // parent 1 heterozygous at allele 99
    parents.getSequenceSpace().flip( 3, 77 );   // parent 1 homozygous at allele 77
    parents.getSequenceSpace().flip( 3, 88 );   // parent 1 heterozygous at allele 88

    xover.update( &parents, mates, &offspring );

    bool obs_state = offspring.getSequenceSpace()(0, 28);        // offspring should be set at 28
    BOOST_REQUIRE_MESSAGE( obs_state, "Unexpected state for offspring at allele 28");

    obs_state = offspring.getSequenceSpace()(1, 77);
    BOOST_REQUIRE_MESSAGE( obs_state, "Unexpected state for offspring at allele 77");

    obs_state = offspring.getSequenceSpace()(1, 33 );
    BOOST_REQUIRE_MESSAGE( !obs_state, "Unexpected state for offspring at allele 33" );

    obs_state = offspring.getSequenceSpace()( 1, 88 );  // after recombination we expect this to be set
    BOOST_REQUIRE_MESSAGE( obs_state, "Unexpected state for offspring at allele 88" );

    obs_state = offspring.getSequenceSpace()( 1, 99 );  // after recombination we expect this to be unset
    BOOST_REQUIRE_MESSAGE( !obs_state, "Unexpected state for offspring at allele 99" );
}

BOOST_AUTO_TEST_SUITE_END()
