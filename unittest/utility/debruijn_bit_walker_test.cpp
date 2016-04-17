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

#include "clotho/utility/debruijn_bit_walker.hpp"
#include "clotho/utility/popcount.hpp"
#include <vector>

template < class Type >
void validate_set_bits( std::vector< unsigned int > & indices, Type origin ) {

    std::vector< unsigned int >::iterator it = indices.begin();

    while( it != indices.end() ) {
        Type obs_bit = ((Type)1 << *it);

        BOOST_REQUIRE_MESSAGE( (obs_bit & origin), "Unexpected set bit index: Observed bit: " << obs_bit );

        origin ^= obs_bit;
        ++it;
    }

    BOOST_REQUIRE_MESSAGE( origin == 0, "Unexpected set bits remaining: Observed: " << origin << "; Expected: 0" );
}

BOOST_AUTO_TEST_SUITE( test_utility )

BOOST_AUTO_TEST_CASE( test_debruijn_walker_char ) {
    typedef char block_type;
    typedef clotho::utility::debruijn_bit_walker< block_type > walker_type;

    block_type origin = 0x13;

    block_type b = origin;
    size_t exp_size = popcount( b );

    std::vector< unsigned int > indices;

    while( b ) {
        unsigned int idx = walker_type::unset_next_index( b );
        indices.push_back( idx );
    }

    size_t obs_size = indices.size();

    BOOST_REQUIRE_MESSAGE( obs_size == exp_size, "Unexpected number of set bits returned; Observed: " << obs_size << "; Expected: " << exp_size );

    validate_set_bits( indices, origin );
}


BOOST_AUTO_TEST_CASE( test_debruijn_walker_uchar ) {
    typedef unsigned char block_type;
    typedef clotho::utility::debruijn_bit_walker< block_type > walker_type;

    block_type origin = 0x13;

    block_type b = origin;
    size_t exp_size = popcount( b );

    std::vector< unsigned int > indices;

    while( b ) {
        unsigned int idx = walker_type::unset_next_index( b );
        indices.push_back( idx );
    }

    size_t obs_size = indices.size();

    BOOST_REQUIRE_MESSAGE( obs_size == exp_size, "Unexpected number of set bits returned; Observed: " << obs_size << "; Expected: " << exp_size );

    validate_set_bits( indices, origin );
}

BOOST_AUTO_TEST_CASE( test_debruijn_walker_unlonglong ) {
    typedef unsigned long long block_type;
    typedef clotho::utility::debruijn_bit_walker< block_type > walker_type;

    block_type origin = 0x1000F33210220307;

    block_type b = origin;
    size_t exp_size = popcount( b );

    std::vector< unsigned int > indices;

    while( b ) {
        unsigned int idx = walker_type::unset_next_index( b );
        indices.push_back( idx );
    }

    size_t obs_size = indices.size();

    BOOST_REQUIRE_MESSAGE( obs_size == exp_size, "Unexpected number of set bits returned; Observed: " << obs_size << "; Expected: " << exp_size );

    validate_set_bits( indices, origin );
}

BOOST_AUTO_TEST_SUITE_END()
