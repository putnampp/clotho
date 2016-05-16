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
#include "clotho/data_spaces/selection/selection.hpp"

#include <boost/random/mersenne_twister.hpp>
#include <boost/property_tree/ptree.hpp>

using namespace clotho::genetics;

typedef double  position_type;
typedef double  weight_type;
typedef unsigned long long block_type;

typedef boost::random::mt19937                                          random_engine_type;
typedef qtl_allele_vectorized< position_type, weight_type >             allele_type;
typedef genetic_space< allele_type, block_type, clotho::genetics::column_aligned >                        genetic_space_type;

BOOST_AUTO_TEST_SUITE( test_data_space )

BOOST_AUTO_TEST_CASE( selection_generator_test ) {
    typedef SelectionGenerator< random_engine_type, genetic_space_type >    selection_generator_type;
    random_engine_type  rand;
    boost::property_tree::ptree config;

    selection_generator_type sel( &rand, config );
}

BOOST_AUTO_TEST_CASE( fitness_selection_generator_test ) {
    typedef SelectionGenerator< random_engine_type, fitness_selection< genetic_space_type > > selection_generator_type;
    random_engine_type  rand;
    boost::property_tree::ptree config;

    selection_generator_type sel( &rand, config );
}

BOOST_AUTO_TEST_CASE( assortative_selection_generator_test ) {
    typedef SelectionGenerator< random_engine_type, assortative_selection< genetic_space_type > > selection_generator_type;
    random_engine_type  rand;
    boost::property_tree::ptree config;

    selection_generator_type sel( &rand, config );
}

BOOST_AUTO_TEST_SUITE_END()
