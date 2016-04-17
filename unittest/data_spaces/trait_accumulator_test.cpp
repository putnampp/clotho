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

#include "clotho/data_spaces/allele_space/qtl_allele.hpp"
#include "clotho/data_spaces/population_space/genetic_space.hpp"
#include "clotho/data_spaces/phenotype_evaluator/trait_accumulator.hpp"

BOOST_AUTO_TEST_SUITE( test_data_space )

BOOST_AUTO_TEST_CASE( trait_accum_size_test ) {

    typedef clotho::genetics::qtl_allele_vectorized< double, double >   allele_type;
    typedef unsigned long long                                          block_type;
    typedef clotho::genetics::genetic_space< allele_type, block_type >  genetic_space_type;
    typedef clotho::genetics::trait_accumulator< genetic_space_type >   trait_accumulator_type;

    size_t exp_alleles = 12, exp_genomes = 110, exp_traits = 15;
    genetic_space_type  gs;

    gs.grow( exp_genomes, exp_alleles );
    gs.getAlleleSpace().grow( exp_alleles, exp_traits );
    
    trait_accumulator_type weigths( gs );
}

BOOST_AUTO_TEST_SUITE_END()
