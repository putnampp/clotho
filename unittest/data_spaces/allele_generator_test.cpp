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

BOOST_AUTO_TEST_CASE( base_allele_generator_test ) {
    typedef base_allele_vectorized< position_type > allele_type;
    typedef AlleleGenerator< rng_type, allele_type > allele_generator_type;

    rng_type    rng;
    boost::property_tree::ptree exp_config;
    boost::property_tree::ptree obs_config;

    allele_generator_type   all_gen( &rng, obs_config );

    BOOST_REQUIRE_MESSAGE( obs_config == exp_config, "Unexpected configuration of default base allele" );
}

BOOST_AUTO_TEST_CASE( neutral_allele_generator_test ) {
    typedef neutral_allele_vectorized< position_type > allele_type;
    typedef AlleleGenerator< rng_type, allele_type > allele_generator_type;

    rng_type    rng;
    boost::property_tree::ptree exp_config;

    exp_config.put( "neutral.p", 0.5 );     // qtl alleles may be neutral

    boost::property_tree::ptree obs_config;

    allele_generator_type   all_gen( &rng, obs_config );

    BOOST_REQUIRE_MESSAGE( obs_config == exp_config, "Unexpected configuration of default neutral allele" );
}

BOOST_AUTO_TEST_CASE( qtl_allele_generator_test ) {
    typedef qtl_allele_vectorized< position_type, weight_type > allele_type;
    typedef AlleleGenerator< rng_type, allele_type > allele_generator_type;

    rng_type    rng;
    boost::property_tree::ptree exp_config;

    exp_config.put( "neutral.p", 0.5 );     // qtl alleles may be neutral
    exp_config.put( "traits.mean", 0.0 );
    exp_config.put( "traits.sigma", 1.0 );

    boost::property_tree::ptree obs_config;

    allele_generator_type   all_gen( &rng, obs_config );

    BOOST_REQUIRE_MESSAGE( obs_config == exp_config, "Unexpected configuration default qtl allele" );
}

BOOST_AUTO_TEST_CASE( qtl_allele_generator_test2 ) {
    typedef qtl_allele_vectorized< position_type, weight_type > allele_type;
    typedef AlleleGenerator< rng_type, allele_type > allele_generator_type;

    rng_type    rng;
    boost::property_tree::ptree config;

    allele_type all( 2, 3);

    allele_generator_type   all_gen( &rng, config );

    all_gen.generate( all, 1 );
}

BOOST_AUTO_TEST_CASE( qtl_allele_push_back_test ) {
    typedef qtl_allele_vectorized< position_type, weight_type > allele_type;
    typedef AlleleGenerator< rng_type, allele_type >            allele_generator_type;

    typedef typename allele_type::trait_iterator                trait_iterator;

    rng_type    rng;
    boost::property_tree::ptree config;

    allele_generator_type agen( &rng, config );

    size_t exp_alleles = 350, exp_traits = 10;
    allele_type gs, fixed;

    gs.grow( exp_alleles, exp_traits );

    size_t i = 0;
    while( i < exp_alleles ) {
        agen( gs, i++);
    }

    size_t exp_fixed_idx = 220;   

    fixed.push_back( gs, exp_fixed_idx );

    size_t obs_fixed_count = fixed.allele_count();

    BOOST_REQUIRE_MESSAGE( obs_fixed_count == 1, "Unexpected number of fixed alleles; Observed: " << obs_fixed_count << "; Expected: 1" );

    position_type exp_pos = gs.getPositionAt( exp_fixed_idx );
    position_type obs_pos = fixed.getPositionAt( 0 );

    BOOST_REQUIRE_MESSAGE( exp_pos == obs_pos, "Unexpected fixed allele position; Observed: " << obs_pos << "; Expected: " << exp_pos );

    bool exp_neutral = gs.getNeutralAt( exp_fixed_idx );
    bool obs_neutral = fixed.getNeutralAt( 0 );

    BOOST_REQUIRE_MESSAGE( exp_neutral == obs_neutral, "Unexpected fixed allele neutrality; Observed: " << obs_neutral << "; Expected: " << exp_neutral );

    size_t obs_traits = fixed.trait_count();

    BOOST_REQUIRE_MESSAGE( exp_traits == obs_traits, "Unexpected fixed allele trait count; Observed: " << obs_traits << "; Expected: " << exp_traits );

    trait_iterator exp_it = gs.getTraitIterator( exp_fixed_idx );
    trait_iterator obs_it = fixed.getTraitIterator( 0 );

    while( obs_it.hasNext() && exp_it.hasNext() ) {
        weight_type exp = exp_it.next();
        weight_type obs = obs_it.next();

        BOOST_REQUIRE_MESSAGE( exp == obs, "Unexpected fixed allele trait weight; Observed: " << obs << "; Expected: " << exp );
    }   
}

BOOST_AUTO_TEST_SUITE_END()
