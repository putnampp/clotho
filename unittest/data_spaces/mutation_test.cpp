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
#include "clotho/data_spaces/mutation/mutation_generator.hpp"

#include "clotho/data_spaces/free_space/free_space.hpp"

#include <boost/property_tree/ptree.hpp>
#include <boost/random/mersenne_twister.hpp>

typedef boost::random::mt19937  random_engine_type;

typedef double                                                                          position_type;
typedef double                                                                          weight_type;

typedef unsigned long long                                                              block_type;

typedef clotho::genetics::qtl_allele_vectorized< position_type, weight_type >           allele_type;
typedef clotho::genetics::genetic_space< allele_type, block_type, clotho::genetics::column_aligned > genetic_space_type;
typedef clotho::genetics::MutationGenerator< random_engine_type, genetic_space_type >   mutation_generator_type;

typedef clotho::genetics::FreeSpaceAnalyzer< genetic_space_type >                       free_space_type;

BOOST_AUTO_TEST_SUITE( test_data_space )

BOOST_AUTO_TEST_CASE( qtl_allele_mutator_test ) {

    boost::property_tree::ptree config;
    random_engine_type  rand(1234);

    size_t  exp_genomes = 20, exp_alleles = 200;

    genetic_space_type  gs;
    mutation_generator_type mutate( &rand, config );

    gs.grow( exp_genomes, exp_alleles );

    gs.getSequenceSpace().clear();

    size_t obs_free_size = gs.getAlleleSpace().free_size();

    BOOST_REQUIRE_MESSAGE( obs_free_size == exp_alleles, "Unexpected number of free alleles within the population; Observed: " << obs_free_size << "; Expected: " << exp_alleles );

    size_t  exp_mutations = 15;

    mutate( &gs, exp_mutations );

    // check that free space was reduced as expected
    obs_free_size = gs.getAlleleSpace().free_size();

    BOOST_REQUIRE_MESSAGE( obs_free_size == (exp_alleles - exp_mutations), "Unexpected number of free alleles within the population; Observed: " << obs_free_size << "; Expected: " << (exp_alleles - exp_mutations) );

    free_space_type fsa( gs );

    obs_free_size = fsa.free_size();

    BOOST_REQUIRE_MESSAGE( obs_free_size == (exp_alleles - exp_mutations), "Unexpected number of free alleles analyzed; Observed: " << obs_free_size << "; Expected: " << (exp_alleles - exp_mutations) );

}

BOOST_AUTO_TEST_SUITE_END()
