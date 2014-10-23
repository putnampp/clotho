#include <boost/test/unit_test.hpp>

#include <vector>
#include <type_traits>

#include "clotho/utility/bit_walker.hpp"

typedef unsigned long block_type;
typedef unsigned short sub_block_type;
typedef clotho::utility::block_walker< block_type, sub_block_type > walker_type;

typedef clotho::utility::block_walker< block_type, sub_block_type, clotho::utility::bit_block_walker< sizeof(block_type) * 8> > other_walker_type;

struct set_bit_vector_op {

    set_bit_vector_op() {}

    void operator()( unsigned int idx ) {
        indices.push_back(idx);
    }

    void reset() { indices.clear(); }

    std::vector< unsigned int > indices;
};

BOOST_AUTO_TEST_SUITE( test_utility )

/**
 * Verify bits_per_block are initialized appropriately
 */
BOOST_AUTO_TEST_CASE( test_bit_walker_bit_sizes ) {

    const unsigned int bpb = walker_type::bits_per_block;
    const unsigned int bpsb = walker_type::bits_per_subblock;
    BOOST_REQUIRE_MESSAGE( bpb == sizeof(block_type) * 8, "Unexpected bits_per_block: " << bpb << " (" << (sizeof(block_type) * 8) << ")");
    BOOST_REQUIRE_MESSAGE( bpsb == sizeof(sub_block_type) * 8, "Unexpected bits_per_subblock: " << bpsb << " (" << (sizeof(sub_block_type) * 8) << ")");
}

BOOST_AUTO_TEST_CASE( test_bit_node_init ) {
    walker_type::bit_walker_type bw;

    for( unsigned int i = 1; i < walker_type::bit_walker_type::max_nodes; ++i ) {
        unsigned int s = (i & (1 << bw.getNode()[i].bit_index));
        BOOST_REQUIRE_MESSAGE( s != 0, "Unexpected bit state: " << i << " @ " << bw.getNode()[i].bit_index );

        if( bw.getNode()[i].next ) {
            unsigned int k = (i >> bw.getNode()[i].bit_shift_next);
            BOOST_REQUIRE_MESSAGE( k == bw.getNode()[i].next, "Unexpected next node: " << i << " -> " << k << " (" << bw.getNode()[i].next << ")");
        }
    }
}

/**
 * Initialize a block with every odd bit set, and initialize an expected index vector
 *
 * Apply the set_bit_vector_op
 *
 * Verify that the result vector matches the expected vector
 */
BOOST_AUTO_TEST_CASE( test_bit_walker_set_bits ) {
    std::vector< unsigned int > expected;

    block_type _bits = 0, mask = 1;
    for( unsigned int i = 0; i < walker_type::bits_per_block; ++i ) {
        if( i & 1 ) {
            expected.push_back( i );

            _bits |= mask;
        }
        mask <<= 1;
    }

    set_bit_vector_op sbv;
    walker_type::apply( _bits, sbv );

    BOOST_REQUIRE_MESSAGE( sbv.indices == expected, "Unexpected list of set bits was returned" );
}

/**
 * Initialize a block with every odd bit set, and initialize an expected index vector
 *
 * Apply the set_bit_vector_op
 *
 * Verify that the result vector matches the expected vector
 */
BOOST_AUTO_TEST_CASE( test_other_bit_walker_set_bits ) {
    BOOST_REQUIRE_MESSAGE( (std::is_same< other_walker_type::bit_walker_type, clotho::utility::bit_block_walker< sizeof(block_type) * 8 > >::value), "Unexpected bit block walker (test 1)" );
    BOOST_REQUIRE_MESSAGE( (!std::is_same< other_walker_type::bit_walker_type, clotho::utility::bit_block_walker< 16 > >::value), "Unexpected bit block walker (test 2)" );

    std::vector< unsigned int > expected;

    block_type _bits = 0, mask = 1;
    for( unsigned int i = 0; i < walker_type::bits_per_block; ++i ) {
        if( i & 1 ) {
            expected.push_back( i );

            _bits |= mask;
        }
        mask <<= 1;
    }

    set_bit_vector_op sbv;
    other_walker_type::apply( _bits, sbv );

    BOOST_REQUIRE_MESSAGE( sbv.indices == expected, "Unexpected list of set bits was returned" );
}

BOOST_AUTO_TEST_SUITE_END()
