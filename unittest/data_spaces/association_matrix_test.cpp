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

#include "clotho/data_spaces/association_matrix.hpp"
#include "clotho/utility/bit_helper.hpp"

using namespace clotho::genetics;

BOOST_AUTO_TEST_SUITE( test_data_space )

BOOST_AUTO_TEST_CASE( association_matrix_sizing_test ) {

    typedef unsigned long long block_type;
    typedef clotho::utility::BitHelper< block_type > bit_helper_type;

    size_t exp_rows = 5, exp_cols = 2;

    size_t exp_size = exp_rows * exp_cols;
    size_t exp_alloc_size = (exp_cols / bit_helper_type::BITS_PER_BLOCK + 1) * exp_rows;

    association_matrix< block_type > amat( exp_rows, exp_cols );

    size_t obs_size = amat.size();
    size_t obs_alloc_size = amat.allocated_size();
    
    BOOST_REQUIRE_MESSAGE( obs_size == exp_size, "Unexpected matrix size; Observed size: " << obs_size << "; Expected size: " << exp_size );
    BOOST_REQUIRE_MESSAGE( obs_alloc_size == exp_alloc_size, "Unexpected allocated size: Observed allocated size: " << obs_alloc_size << "; Expected allocated size: " << exp_alloc_size );

}

BOOST_AUTO_TEST_CASE( association_matrix_flip_test ) {
    typedef unsigned long long block_type;
    typedef clotho::utility::BitHelper< block_type > bit_helper_type;

    size_t exp_rows = 5, exp_cols = 12;

    association_matrix< block_type > amat( exp_rows, exp_cols );

    size_t flip_row = 3, flip_col = 6;

    bool exp_state = !amat( flip_row, flip_col );

    amat.flip(flip_row, flip_col );
    bool obs_state = amat( flip_row, flip_col );

    BOOST_REQUIRE_MESSAGE( obs_state == exp_state, "Unexpected flip state; Observed state: " << obs_state << "; Expected state: " << exp_state );
}

BOOST_AUTO_TEST_SUITE_END()
