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
#include "clotho/data_spaces/free_space/free_space.hpp"

typedef double              position_type;
typedef double              weight_type;
typedef unsigned long long  block_type;
typedef clotho::genetics::qtl_allele_vectorized< position_type, weight_type >       allele_type;
typedef clotho::genetics::genetic_space< allele_type, block_type >                  genetic_space_type;
typedef clotho::genetics::FreeSpaceAnalyzer< genetic_space_type >                   free_space_type;

BOOST_AUTO_TEST_SUITE( test_data_space )

BOOST_AUTO_TEST_CASE( free_space_test ) {
    size_t  exp_alleles = 200, exp_genomes = 11;
    genetic_space_type  gs;

    gs.grow( exp_genomes, exp_alleles );

    gs.getSequenceSpace().clear();
    
    // fix index
    size_t fixed = 101;
    size_t i = 0;
    while( i < 2 * exp_genomes) {
        gs.getSequenceSpace().flip( i, fixed );
        ++i;
    }

    // make random allele (33) variable
    gs.getSequenceSpace().flip( 2, 33 );

    size_t exp_fixed = 1, exp_var = 1, exp_free = exp_alleles - exp_var, exp_lost = exp_free - exp_fixed;

    free_space_type fs;

    fs.update( gs );

    size_t obs_fixed = fs.fixed_size();

    BOOST_REQUIRE_MESSAGE( obs_fixed == exp_fixed, "Unexpected number of fixed alleles; Observed: " << obs_fixed << "; Expected: " << exp_fixed );

    size_t obs_lost = fs.lost_size();

    BOOST_REQUIRE_MESSAGE( obs_lost == exp_lost, "Unexpected number of lost alleles; Observed: " << obs_lost << "; Expected: " << exp_lost );

    size_t obs_free = fs.free_size();

    BOOST_REQUIRE_MESSAGE( obs_free == exp_free, "Unexpected number of free alleles; Observed: " << obs_free << "; Expected: " << exp_free );
}

BOOST_AUTO_TEST_SUITE_END()
