#include <boost/test/unit_test.hpp>

#include "test_element.h"
#include "clotho/powerset/variable_subset.hpp"

typedef unsigned long Block;

typedef clotho::powersets::block_map< test_element, Block, clotho::powersets::normalized_key< test_element > > bmap;
typedef clotho::powersets::variable_subset< test_element, Block, bmap > subset_type;
typedef typename subset_type::powerset_type powerset_type;

BOOST_AUTO_TEST_SUITE( test_powerset )

/**
 *  Test that the expected width per block is initialized correctly
 *
 */
BOOST_AUTO_TEST_CASE( test_bmap ) {
    double bwidth = bmap::width_per_bin;
    BOOST_REQUIRE_MESSAGE( bwidth == (1.0/64.0), "Unexpected width per bin: " << bwidth << "(" << (1.0/64.0) << ")");
}

/**
 * Create an empty powerset
 *
 * Manually append element to powerset
 * Verify that new element is assigned to expected index
 *
 * Verify allocated size
 */
BOOST_AUTO_TEST_CASE( create_powerset_append) {
    powerset_type ps;

    typename powerset_type::element_index_type idx = ps.appendElement( test_element( 0., 1.0) );

    BOOST_REQUIRE_MESSAGE( idx == 0, "Unexpected index " << idx << " returned for " << 0 << "(" << 0 << ")" );
    BOOST_REQUIRE_MESSAGE( ps.variable_allocated_size() == bmap::bits_per_block, "Unexpected variable space: " << ps.variable_allocated_size() );

    ps.updateFreeIndex( idx, false );
    BOOST_REQUIRE_MESSAGE( ps.free_size() == (bmap::bits_per_block - 1), "Unexpected number of free space: " << ps.free_size() );
}

/**
 * Create an empty powerset
 *
 * Manually add element to powerset (results in append)
 * Verify that new element is assigned to expected index
 *
 * Verify allocated size
 */
BOOST_AUTO_TEST_CASE( create_powerset_add) {

    powerset_type ps;

    typename powerset_type::element_index_type idx = ps.addElement( test_element( 0., 1.0) );
    BOOST_REQUIRE_MESSAGE( idx == 0, "Unexpected index " << idx << " returned for " << 0 << "(" << 0 << ")" );

    BOOST_REQUIRE_MESSAGE( ps.size() == 1, "Unexpected size" );
    BOOST_REQUIRE_MESSAGE( ps.variable_allocated_size() == bmap::bits_per_block, "Unexpected variable space: " << ps.variable_allocated_size() );
}

/**
 * Create an empty powerset
 *
 * Manually append first element to powerset
 * Verify that new element is assigned to expected index
 *
 * add the remaining elements of a block
 * Verify that each new element is assigned to expected index
 *
 * Verify allocated size
 */
BOOST_AUTO_TEST_CASE( create_powerset_width ) {
    powerset_type ps;

    BOOST_REQUIRE_MESSAGE( ps.empty(), "Unexpected size" );

    typename powerset_type::element_index_type idx = ps.appendElement( test_element( 0., 1.0) );
    
    BOOST_REQUIRE_MESSAGE( idx == 0, "Unexpected index " << idx << " returned for " << 0 << "(" << 0 << ")" );
    BOOST_REQUIRE_MESSAGE( ps.variable_allocated_size() == bmap::bits_per_block, "Unexpected variable space: " << ps.variable_allocated_size() );

    ps.updateFreeIndex( idx, false );
    BOOST_REQUIRE_MESSAGE( ps.free_size() == (bmap::bits_per_block - 1), "Unexpected number of free space: " << ps.free_size() );

    for( unsigned int i = 1; i < bmap::bits_per_block; ++i ) {
        double v = (double) i;
        double k = v / (double) bmap::bits_per_block;
        test_element te(k, 1.0);

        typename powerset_type::element_index_type e_idx = ps.findFreeIndex(te);

        BOOST_REQUIRE_MESSAGE( e_idx == i, e_idx << " != " << i );

        typename powerset_type::element_index_type idx = ps.addElement(te);

        BOOST_REQUIRE_MESSAGE( idx == i, "Unexpected index " << idx << " returned for " << i << "(" << k << "; " << e_idx << ")" );
    }

    const unsigned int bpb = bmap::bits_per_block;
    BOOST_REQUIRE_MESSAGE( ps.variable_allocated_size() == bpb, "Unexpected size: " << ps.variable_allocated_size() << "(" << bpb << ")" );
}

/**
 * Create an empty powerset
 *
 * find_or_create elements in a round robin order.
 * Verify that each new element is assigned to expected index
 */
BOOST_AUTO_TEST_CASE( create_powerset_find_or_create_unique ) {
    powerset_type ps;

    const unsigned int sub_div = 3;
    const double div_offset = (1.0/ ((double)sub_div*bmap::bits_per_block));

    for( unsigned int i = 0; i < sub_div * bmap::bits_per_block; ++i ) {
        double v = (double) (i % bmap::bits_per_block);

        unsigned int s = (i / bmap::bits_per_block );
        
        double k = v / (double) bmap::bits_per_block + (double)s*div_offset;
        test_element te(k, 1.0);

        typename powerset_type::element_index_type idx = ps.find_or_create(te);

        BOOST_REQUIRE_MESSAGE( idx == i, "Unexpected index " << idx << " returned for " << i << "(" << k << ")" );
    }

    const unsigned int bpb = sub_div * bmap::bits_per_block;
    BOOST_REQUIRE_MESSAGE( ps.variable_allocated_size() == bpb, "Unexpected size: " << ps.variable_allocated_size() << "(" << bpb << ")" );
}

/**
 * Create an empty powerset
 *
 * find_or_create an element to the powerset (results in create)
 *
 * Verify that it is in expected position
 *
 * find_or_create different element with same key (results in find)
 *
 * Verify that element index returned is same as earlier element.
 */
BOOST_AUTO_TEST_CASE( create_powerset_find_or_create_collision ) {
    powerset_type ps;

    double k = 0.5;
    test_element te(k, 1.0);
    typename powerset_type::element_index_type idx = ps.find_or_create(te);
    BOOST_REQUIRE_MESSAGE( idx == bmap()(te), "Unexpected index " << idx << " returned for " << 0 << "(" << k << ")" );

    test_element te2(k, 2.0);
    typename powerset_type::element_index_type cidx = ps.find_or_create(te2);
    BOOST_REQUIRE_MESSAGE( cidx == idx, "Unexpected index " << idx << " returned for " << cidx << "(" << k << ")" );
}

/**
 * Create an empty powerset
 *
 * find_or_create elements in a bit position order.
 * Every element after the first in each sub-division should result in a collision.
 *
 * Verify that each new element is assigned to expected index
 */
BOOST_AUTO_TEST_CASE( create_powerset_find_or_create_order_collision ) {
    powerset_type ps;

    const unsigned int sub_div = 3;
    const double div_offset = (1.0/ ((double)sub_div*bmap::bits_per_block));

    for( unsigned int i = 0; i < bmap::bits_per_block; ++i ) {
        for( unsigned int j = 0; j < sub_div; ++j ) {
            double v = (double)i;
            double k = v / (double) bmap::bits_per_block + (double)j*div_offset;

            test_element te(k, 1.0);

            unsigned int e_idx = i + j * bmap::bits_per_block;
            typename powerset_type::element_index_type idx = ps.find_or_create(te);

            BOOST_REQUIRE_MESSAGE( idx == e_idx, "Unexpected index " << idx << " returned for " << e_idx << "(" << k << ")" );
        }
    }

    const unsigned int e_size = sub_div * bmap::bits_per_block;
    BOOST_REQUIRE_MESSAGE( ps.variable_allocated_size() == e_size, "Unexpected size: " << ps.variable_allocated_size() << "(" << e_size << ")" );
}

/**
 * Creates a subset of the powerset
 *
 * Verify that subset was added to family
 *
 * Adds a element
 *
 * Checks that element has been added to the set and subset
 * Checks subset count
 * Verifies expected sizes
 * Verifies expected bit state
 *
 * Removes earlier element
 * Checks subset count 
 * Checks element bit state
 */
BOOST_AUTO_TEST_CASE( create_powerset_subset ) {
    powerset_type ps;

    typename powerset_type::subset_ptr c = ps.create_subset();

    BOOST_REQUIRE_MESSAGE( ps.family_size() == 1, "Subset was not added to family");

    double k = 0.5;
    test_element te(k, 1.0);

    c->addElement( te );

    size_t cnt = c->count();
    BOOST_REQUIRE_MESSAGE( cnt == 1, "Element has not been added to the subset: " << cnt << "(1)" );

    const unsigned int bpb = bmap::bits_per_block;
    BOOST_REQUIRE_MESSAGE( ps.variable_allocated_size() == bpb, "Unexpected size: " << ps.variable_allocated_size() << "(" << bpb << ")" );

    typename powerset_type::element_index_type idx = ps.find_or_create(te);
    BOOST_REQUIRE_MESSAGE( idx == bmap()(te), "Unexpected index " << idx << " returned for " << 0 << "(" << k << ")" );

    BOOST_REQUIRE_MESSAGE( c->check_state( idx ), "Unexpected state for element after add" );

    c->removeElement( te );

    cnt = c->count();
    BOOST_REQUIRE_MESSAGE( cnt == 0, "Element has not been removed to the subset: " << cnt << "(0)" );
    BOOST_REQUIRE_MESSAGE( !c->check_state( idx ), "Unexpected state for element after remove" );
}

/**
 * Create a powerset
 *
 * Create a subset with one element in it, in a sub-scope;
 * Assuming subset_ptr is a shared_ptr, the change in scope should dereference
 * a single instance to the subset, leaving a single instance of the subset
 * in the family.  If the family contains only a single instance of a subset,
 * prune the powerset should remove the subset, and free up lost|fixed
 * elements.
 *
 */
BOOST_AUTO_TEST_CASE( create_powerset_prunespace ) {
    powerset_type ps;

    // change in scope should release one instance of subset_ptr
    {
        typename powerset_type::subset_ptr c = ps.create_subset();

        double k = 0.5;
        test_element te(k, 1.0);

        c->addElement( te );
    }

    BOOST_REQUIRE_MESSAGE( ps.family_size() == 1, "Subset was not added to family");

    const unsigned int bpb = bmap::bits_per_block;
    size_t fs = ps.free_size();
    BOOST_REQUIRE_MESSAGE( fs == (bpb - 1), "Unexpected free size " << fs << "(" << (bpb - 1) << ")" );
    BOOST_REQUIRE_MESSAGE( ps.variable_allocated_size() == bpb, "Unexpected size: " << ps.variable_allocated_size() << "(" << bpb << ")" );

    ps.pruneSpace();

    BOOST_REQUIRE_MESSAGE( ps.family_size() == 0, "Subset was not removed from family");

    fs = ps.free_size();
    BOOST_REQUIRE_MESSAGE( fs == bpb, "Unexpected free size " << fs << "(" << bpb << ")" );
}

/**
 * Create a powerset
 *
 * Create a subset with one element in it, in a sub-scope;
 * Create a second subset with the same element (a duplicate of the first subset)
 *
 * As before the first subset should be dereferenced, and will be removed from
 * the family.  The second subset, however, will not because two references to the
 * subset should exist (one in this scope, and one in the family).
 *
 * Therefore, pruning the powerset space should result in the subsets being 
 * reduced to 1.  Furthermore, since there is a single subset, with a single
 * element, that element becomes fixed within the powerset.  The fixed element
 * should be copied to the fixed subset, and its positions should be freed
 *
 */
BOOST_AUTO_TEST_CASE( create_powerset_prunespace2 ) {
    powerset_type ps;

    double k = 0.5;
    test_element te(k, 1.0);

    // change in scope should release one instance of subset_ptr
    {
        typename powerset_type::subset_ptr c = ps.create_subset();

        c->addElement( te );
    }

    typename powerset_type::subset_ptr c2 = ps.create_subset();
    c2->addElement(te);

    BOOST_REQUIRE_MESSAGE( ps.family_size() == 2, "Subset was not added to family");

    typename powerset_type::element_index_type idx = ps.find_or_create(te);

    const unsigned int bpb = bmap::bits_per_block;
    size_t fs = ps.free_size();
    BOOST_REQUIRE_MESSAGE( fs == (bpb - 1), "Unexpected free size " << fs << "(" << (bpb - 1) << ")" );
    BOOST_REQUIRE_MESSAGE( ps.variable_allocated_size() == bpb, "Unexpected size: " << ps.variable_allocated_size() << " (" << bpb << ")" );

    ps.pruneSpace();

    fs = ps.free_size();
    BOOST_REQUIRE_MESSAGE( ps.family_size() == 1, "Expected subset was not removed from family " << ps.family_size() << " (1)" );
    BOOST_REQUIRE_MESSAGE( fs == bpb, "Unexpected free size " << fs << " (" << bpb << ")" );
    BOOST_REQUIRE_MESSAGE( ps.fixed_size() == 1, "Unexpected fixed size " << ps.fixed_size() << " (1)"); 

    typename powerset_type::element_index_type idx2 = ps.find_or_create(te);
    BOOST_REQUIRE_MESSAGE( idx != idx2, "After pruning, fixed element has same index.  Should be different");

    BOOST_REQUIRE_MESSAGE( idx2 == ps.encode_index(0, true), "Unxpected fixed index " << idx2 << " (" << ps.encode_index(0, true) << ")");

    bool is_fixed  = ps.decode_index(idx2);
    BOOST_REQUIRE_MESSAGE( is_fixed, "Unexpected state decoding false (true)" );
    BOOST_REQUIRE_MESSAGE( idx2 == 0, "Unexpected index decoded " << idx2 << " (0)");
}
BOOST_AUTO_TEST_SUITE_END()
