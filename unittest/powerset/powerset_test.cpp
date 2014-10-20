#include <boost/test/unit_test.hpp>

//#include "clotho/powerset/powerset.hpp"
//
#include "clotho/powerset/block_map.hpp"
#include "clotho/powerset/variable_subset.hpp"

struct test_element {
    double k, v;

    test_element( double _k = 0., double _v = 0.) : k (_k), v(_v) {}
    test_element( const test_element & t ) : k(t.k), v(t.v) {}

    friend bool operator==( const test_element & lhs, const test_element & rhs ); 
    friend bool operator!=( const test_element & lhs, const test_element & rhs ); 
};

inline bool operator==( const test_element & lhs, const test_element & rhs ) {
    return ( lhs.k == rhs.k && lhs.v == rhs.v);
}

inline bool operator!=( const test_element & lhs, const test_element & rhs ) {
    return ( lhs.k != rhs.k && lhs.v != rhs.v);
}

namespace clotho {
namespace powersets {

template <>
struct element_key_of< test_element > {
    typedef double key_type;

    inline key_type operator()( const test_element & t ) { return t.k; }
};

template< class Block >
struct block_map< test_element, Block > {
    typedef test_element    element_type;
    typedef Block           size_type;

    static const unsigned int bits_per_block = sizeof(size_type) * 8;

    inline size_type operator()( const element_type & elem ) {
        assert( 0. <= elem.k && elem.k < 1. );

        return elem.k * bits_per_block;
    }
};

}   // namespace powersets
}   // namespace clotho

typedef clotho::powersets::variable_subset< test_element > subset_type;
typedef typename subset_type::powerset_type powerset_type;
typedef powerset_type::block_map_type bmap;

BOOST_AUTO_TEST_SUITE( test_powerset )

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

BOOST_AUTO_TEST_CASE( create_powerset_prunespace ) {
// TODO
}

BOOST_AUTO_TEST_SUITE_END()
