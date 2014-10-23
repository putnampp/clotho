#include <boost/test/unit_test.hpp>

#include "test_element.h"
#include "clotho/powerset/variable_subset.hpp"
#include "clotho/powerset/variable_subset_recombination.hpp"

typedef unsigned long Block;

typedef clotho::powersets::block_map< test_element, Block, clotho::powersets::normalized_key< test_element > > bmap;
typedef clotho::powersets::variable_subset< test_element, Block, bmap > subset_type;
typedef typename subset_type::powerset_type powerset_type;

template < class Block >
struct simple_classifier {
    typedef Block block_type;

    simple_classifier() : m_res(0) {}

    void operator()( unsigned int idx ) {
        if( idx & 1 ) {
            m_res |= (1 << idx);
        }
    }

    void resetResult() {
        m_res = 0;
    }

    block_type getResult() const {
        return m_res;
    }

    block_type m_res;
};

typedef simple_classifier< Block > classifier_type;
typedef clotho::recombine::recombination< subset_type, classifier_type > recombination_type;

BOOST_AUTO_TEST_SUITE( test_recombination )


BOOST_AUTO_TEST_CASE( recombine_create ) {
    powerset_type ps;

    typename powerset_type::subset_ptr p0 = ps.create_subset();
    typename powerset_type::subset_ptr p1 = ps.create_subset();

    const unsigned int sub_div = 3;
    const double div_offset = (1.0/ ((double)sub_div*bmap::bits_per_block));

    for( unsigned int i = 0; i < sub_div * bmap::bits_per_block; ++i ) {
        double v = (double) (i % bmap::bits_per_block);

        unsigned int s = (i / bmap::bits_per_block );
        
        double k = v / (double) bmap::bits_per_block + (double)s*div_offset;
        test_element te(k, 1.0);
        if( i & 1 ) {
            p1->addElement( te );
        } else {
            p0->addElement( te );
        }

        typename powerset_type::element_index_type idx = ps.find(te);

        BOOST_REQUIRE_MESSAGE( idx == i, "Unexpected index " << idx << " returned for " << i << "(" << k << ")" );
    }

    recombination_type rec;
    classifier_type cfier;

    rec( p1, p0, cfier );

    typename powerset_type::subset_ptr c = ps.create_subset( *rec.getResultSequence() );

    BOOST_REQUIRE_MESSAGE( *p1 == *c, "Unexpected result");   
}

BOOST_AUTO_TEST_SUITE_END()
