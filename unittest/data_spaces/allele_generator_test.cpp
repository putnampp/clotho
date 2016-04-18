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

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <iostream>

using namespace clotho::genetics;

typedef boost::random::mt19937  rng_type;
typedef double position_type;
typedef double weight_type;

BOOST_AUTO_TEST_SUITE( test_data_space )

BOOST_AUTO_TEST_CASE( allele_generator_test ) {
    typedef qtl_allele_vectorized< position_type, weight_type > allele_type;
    typedef AlleleGenerator< rng_type, allele_type > allele_generator_type;

    rng_type    rng;
    boost::property_tree::ptree exp_config;

    exp_config.put( "neutral.p", 0.5 );     // qtl alleles may be neutral
    exp_config.put( "traits.mean", 0.0 );
    exp_config.put( "traits.sigma", 1.0 );

    boost::property_tree::ptree obs_config;

    allele_generator_type   all_gen( &rng, obs_config );

    BOOST_REQUIRE_MESSAGE( obs_config == exp_config, "Unexpected configuration" );

    boost::property_tree::write_json( std::cout, obs_config );
    boost::property_tree::write_json( std::cout, exp_config );
}

BOOST_AUTO_TEST_SUITE_END()
