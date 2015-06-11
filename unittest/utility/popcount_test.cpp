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

#include "clotho/utility/popcount.hpp"

#define TEST_PC( x, e) \
    BOOST_REQUIRE_MESSAGE( popcount( x ) == e, "Unexpected popcount: pc(" << x << ") =  pc(0x" << std::hex << x << ") = " << std::dec << popcount(x) << " !=  " << e)

BOOST_AUTO_TEST_SUITE( test_utility )

BOOST_AUTO_TEST_CASE( test_popcount ) {
    int x = -129;   // 0xffffff7f

    TEST_PC( x, 31 );

    char y = -51;   // 0xCD
    TEST_PC( y, 5);

    long z = -129 * 0x0000000100000000;   // 0xffffff7f00000000
    TEST_PC( z, 31);
}

BOOST_AUTO_TEST_SUITE_END()
