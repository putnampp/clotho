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
#include <vector>
#include <iostream>
#include <iomanip>

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

BOOST_AUTO_TEST_CASE( association_matrix_raw_iterator_test ) {
    typedef unsigned long long block_type;

    typedef unsigned long long block_type;
    typedef clotho::utility::BitHelper< block_type >    bit_helper_type;

    typedef association_matrix< block_type >            association_type;
    typedef typename association_type::block_iterator   iterator;

    size_t exp_rows = 5, exp_cols = 240;

    association_type amat( exp_rows, exp_cols );

    size_t exp_blocks = exp_rows * (exp_cols / bit_helper_type::BITS_PER_BLOCK + 1);

    amat.clear(); 

    iterator it = amat.raw_iterator();

    size_t obs_blocks = 0;
    size_t i = 0, j = 0;
    while( it.hasNext() ) {
        block_type b = it.next();

        BOOST_REQUIRE_MESSAGE( b == 0, "Unexpected block state at row: " << i << ", column: " << j << "; Observed: " << b << "; Expected: 0" );

        if( ++j >= it.max_columns() ) {
            j = 0;
            ++i;
        }
        ++obs_blocks;
    }

    BOOST_REQUIRE_MESSAGE( exp_blocks == obs_blocks, "Unexpected number of raw blocks; Observed: " << obs_blocks << "; Expected: " << exp_blocks );
}

BOOST_AUTO_TEST_CASE( association_matrix_raw_iterator_test2 ) {
    typedef unsigned long long block_type;

    typedef unsigned long long block_type;
    typedef clotho::utility::BitHelper< block_type >    bit_helper_type;

    typedef association_matrix< block_type >            association_type;
    typedef typename association_type::block_iterator   iterator;

    size_t exp_rows = 5, exp_cols = 240;

    association_type amat( exp_rows, exp_cols );

    size_t exp_block_count = exp_rows * (exp_cols / bit_helper_type::BITS_PER_BLOCK + 1);

    amat.clear();

    std::vector< block_type > exp_blocks;
    
    size_t i = 0, j = 0, k = 0;
    block_type _b = 0;
    while( i < exp_rows ) {
        amat.flip( i, j );

        if( k != j / bit_helper_type::BITS_PER_BLOCK ) {
            exp_blocks.push_back(_b);
            _b = 0;
            ++k;
        }

        _b |= (1L << (j % bit_helper_type::BITS_PER_BLOCK));

        j += 15;
        if( j >= exp_cols ) {
            exp_blocks.push_back( _b );
            j = 0;
            k = 0;
            _b = 0;
            ++i;
        }
    }

    iterator it = amat.raw_iterator();

    size_t obs_blocks = 0;
    i = 0; j = 0;

    std::vector< block_type >::iterator exp_it = exp_blocks.begin();

    while( it.hasNext() ) {
        block_type b = it.next();

        size_t idx = j * it.max_rows() + i;

        BOOST_REQUIRE_MESSAGE( b == exp_blocks[idx], "Unexpected block state at row: " << std::dec << i << ", column: " << std::dec << j << "; Observed: " << std::hex << std::setw( sizeof( block_type ) * 2) << std::setfill('0') << b << "; Expected [" << idx << "]: " << std::hex << std::setw( sizeof(block_type) * 2) << std::setfill('0') << exp_blocks[idx] );

        if( ++j >= it.max_columns() ) {
            j = 0;
            ++i;
        }
        ++obs_blocks;
    }

    BOOST_REQUIRE_MESSAGE( exp_block_count == obs_blocks, "Unexpected number of raw blocks; Observed: " << obs_blocks << "; Expected: " << exp_block_count );
}

BOOST_AUTO_TEST_SUITE_END()
