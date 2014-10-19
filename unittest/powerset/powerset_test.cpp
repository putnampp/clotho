#include <boost/test/unit_test.hpp>

//#include "clotho/powerset/powerset.hpp"
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

BOOST_AUTO_TEST_CASE( create_powerset_simple) {

    powerset_type ps;

    double k = 0.75;
    test_element te( k, 2.0 );

    typename powerset_type::element_index_type idx = ps.find_or_create(te);
    BOOST_REQUIRE_MESSAGE( idx == bmap()(te), "Unexpected index returned" );
    BOOST_REQUIRE_MESSAGE( ps.size() == 1, "Unexpected size" );
    BOOST_REQUIRE_MESSAGE( ps.variable_allocated_size() == bmap::bits_per_block, "Unexpected variable space: " << ps.variable_allocated_size() );

}

BOOST_AUTO_TEST_CASE( create_powerset_width ) {
    powerset_type ps;

    BOOST_REQUIRE_MESSAGE( ps.empty(), "Unexpected size" );
    for( unsigned int i = 0; i < bmap::bits_per_block; ++i ) {
        double v = (double) i;
        double k = v / (double) bmap::bits_per_block;
        test_element te(k, 1.0);

        typename powerset_type::element_index_type idx = ps.find_or_create(te);

        BOOST_REQUIRE_MESSAGE( idx == i, "Unexpected index " << idx << " returned for " << i << "(" << k << ")" );
    }

    BOOST_REQUIRE_MESSAGE( ps.size() == bmap::bits_per_block, "Unexpected size" );
}

BOOST_AUTO_TEST_SUITE_END()
