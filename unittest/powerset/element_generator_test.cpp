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
