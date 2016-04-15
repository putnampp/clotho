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

#include "clotho/data_spaces/allele_space/neutral_allele.hpp"
#include "clotho/data_spaces/allele_space/qtl_allele.hpp"

using namespace clotho::genetics;

BOOST_AUTO_TEST_SUITE( test_data_space )

BOOST_AUTO_TEST_CASE( base_allele_vectorized_test ) {
    neutral_allele_vectorized< double > all;

    size_t  test_idx = 10;
    bool    test_neutral = true;
    all.setNeutralAt( test_idx, test_neutral );

    size_t s = all.size();

    BOOST_REQUIRE_MESSAGE( s == test_idx + 1, "Unexpected size of vectored alleles; Observed: " << s << "; Expected: " << (test_idx + 1));

    bool obs_neutral = all.getNeutralAt( test_idx );

    BOOST_REQUIRE_MESSAGE( obs_neutral == test_neutral, "Unexpected neutrality at index " << test_idx << "; Observed: " << obs_neutral << "; Expected: " << test_neutral);
}

BOOST_AUTO_TEST_CASE( qtl_allele_vectorized_test ) {

    std::cout << "QTL Allele Test:" << std::endl;

    qtl_allele_vectorized< double, double > all;

    size_t obs_size = all.size();
    size_t obs_alleles = all.allele_count();
    size_t obs_traits = all.trait_count();
    size_t obs_allocated_size = all.allocated_size();

    BOOST_REQUIRE_MESSAGE( obs_size == 1, "Unexpected size upon initialization; Observed: " << obs_size << "; Expected: 1" );
    BOOST_REQUIRE_MESSAGE( obs_alleles == 1, "Unexpected alleles upon initialization; Observed: " << obs_alleles << "; Expected: 1" );
    BOOST_REQUIRE_MESSAGE( obs_traits == 1, "Unexpected traits upon initialization; Observed: " << obs_traits << "; Expected: 1" );
    BOOST_REQUIRE_MESSAGE( obs_allocated_size == 1, "Unexpected allocated upon initialization; Observed: " << obs_allocated_size << "; Expected: 1" );

    size_t exp_alleles = 10;
    all.grow( exp_alleles );

    obs_size = all.size();
    obs_alleles = all.allele_count();
    obs_traits = all.trait_count();
    obs_allocated_size = all.allocated_size();


    BOOST_REQUIRE_MESSAGE( obs_size == exp_alleles, "Unexpected size upon initialization; Observed: " << obs_size << "; Expected: " << exp_alleles );
    BOOST_REQUIRE_MESSAGE( obs_alleles == exp_alleles, "Unexpected alleles upon initialization; Observed: " << obs_alleles << "; Expected: " << exp_alleles );
    BOOST_REQUIRE_MESSAGE( obs_traits == 1, "Unexpected traits upon initialization; Observed: " << obs_traits << "; Expected: 1" );
    BOOST_REQUIRE_MESSAGE( obs_allocated_size == exp_alleles, "Unexpected allocated upon initialization; Observed: " << obs_allocated_size << "; Expected: " << exp_alleles );

    size_t exp_traits = 5;

    all.grow( exp_alleles, exp_traits );

    obs_size = all.size();
    obs_alleles = all.allele_count();
    obs_traits = all.trait_count();
    obs_allocated_size = all.allocated_size();


    BOOST_REQUIRE_MESSAGE( obs_size == exp_alleles, "Unexpected size upon initialization; Observed: " << obs_size << "; Expected: " << exp_alleles );
    BOOST_REQUIRE_MESSAGE( obs_alleles == exp_alleles, "Unexpected alleles upon initialization; Observed: " << obs_alleles << "; Expected: " << exp_alleles );
    BOOST_REQUIRE_MESSAGE( obs_traits == exp_traits, "Unexpected traits upon initialization; Observed: " << obs_traits << "; Expected: " << exp_traits );
    BOOST_REQUIRE_MESSAGE( obs_allocated_size == (exp_alleles * exp_traits), "Unexpected allocated upon initialization; Observed: " << obs_allocated_size << "; Expected: " << (exp_alleles * exp_traits) );

}

BOOST_AUTO_TEST_SUITE_END()
