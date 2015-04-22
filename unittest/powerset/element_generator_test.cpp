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

#include "test_element.h"
#include "clotho/mutation/element_generator.hpp"

typedef clotho::powersets::element_key_of< test_element >::key_type          key_type;

struct test_key_generator {
    key_type operator()( key_type & k ) { return k; }   
};

typedef clotho::mutations::element_generator< test_element, test_key_generator >  elem_generator;

BOOST_AUTO_TEST_SUITE( test_element_generation )

/**
 *  Test that the expected width per block is initialized correctly
 *
 */
BOOST_AUTO_TEST_CASE( key_generator_test ) {
    test_key_generator tkg;
    elem_generator gen( tkg );

    key_type k = 0.25;

    test_element te = gen(k);

    BOOST_REQUIRE_MESSAGE( te.k == k, "Unexpected key: " << te.k << "(" << k << ")" );
}

BOOST_AUTO_TEST_SUITE_END()
