#include <boost/test/unit_test.hpp>

#include "../unittest_config.h"

#include "test_element.h"
#include "clotho/utility/bit_block_iterator.hpp"
#include "clotho/powerset/variable_subset.hpp"
#include "clotho/powerset/variable_subset_recombination.hpp"

#include "clotho/classifiers/region_classifier.hpp"

typedef unsigned long Block;

typedef clotho::powersets::block_map< test_element, Block, clotho::powersets::normalized_key< test_element > > bmap;
typedef clotho::powersets::variable_subset< test_element, Block, bmap > subset_type;
typedef typename subset_type::powerset_type powerset_type;

/**
 * Simple classifier which just classifies all odd indices as being (1 - base) and (0 - alt)
 *
 * The idea being the recombination results as all odd elements from base, and all even
 * elements from alt
 */
struct odd_from_base_classifier {
    odd_from_base_classifier() {}

    template < class ElementIterator >
    bool operator()( ElementIterator it, size_t idx ) const {
        return (idx & 1);
    }

    bool operator()( const test_element & elem ) const {
        return (elem.idx & 1);
    }
};

typedef clotho::recombine::inspection::tag::copy_matching_classify_mismatch inspection_type;
typedef clotho::recombine::walker::tag::inline_dynamic_classify             block_walker_type;
typedef odd_from_base_classifier classifier_type;
typedef clotho::recombine::recombination< subset_type, classifier_type, inspection_type, block_walker_type > recombination_type;

/**
 * A range based classifier.
 *
 * Assumes that recombination_points represent the upper bounds of each key range.
 * For example, the ranges of {0.25, 0.5, 0.75} would then be
 * (-inf, 0.25) -> base
 * [0.25, 0.5) -> alt
 * [0.5, 0.75) -> base
 * [0.75, inf) -> alt
 *
 * Classifier also takes an iterator to a list of elements.  Input index values are
 * relative offsets to the current list element.
 *
 * The result is a bit mask which classifies the list elements at relative indices
 * as being either in a "base" region (1), or "alt" region (0).
 */

typedef clotho::classifiers::region_classifier< test_element > classifier2_type;
typedef clotho::recombine::recombination< subset_type, classifier2_type, inspection_type, block_walker_type > recombination2_type;

BOOST_AUTO_TEST_SUITE( test_recombination )

/**
 * Testing odd_from_base_classifier recombination with empty sequences (NULL shared_ptr)
 *
 * Case 1 - recombines empty sequences, 
 * Case 2 & 3 - recombine empty sequence with sequence with single element (and vice versa).
 * Case 4 - recombine same sequence with element
 */
BOOST_AUTO_TEST_CASE( recombine_empty ) {
    powerset_type ps;

    typename powerset_type::subset_ptr p0 = ps.empty_set();
    typename powerset_type::subset_ptr p1 = ps.empty_set();
    typename powerset_type::subset_ptr c1_exp = ps.create_subset();

    recombination_type rec;
    classifier_type cfier;

    rec( p0, p1, cfier );

    BOOST_REQUIRE_MESSAGE( rec.isEmpty(), "Expected that with p0 and p1 being empty, recombination would be empty");
    BOOST_REQUIRE_MESSAGE( rec.isMatchBase(), "Expected that with p0 and p1 being empty, recombination should match base");
    BOOST_REQUIRE_MESSAGE( rec.isMatchAlt(), "Expected that with p0 and p1 being empty, recombination should match alt");

    unsigned int index = 0;
    test_element te( 0.5, 1.0, index++ );
    c1_exp->addElement( te );

    rec( p0, c1_exp, cfier );

    BOOST_REQUIRE_MESSAGE( !rec.isEmpty(), "Nothing recombined. Expected that with base is empty, even element from c1_exp." );
    BOOST_REQUIRE_MESSAGE( rec.isMatchAlt(), "Expected that with base is empty, even element from c1_exp." );

    rec( c1_exp, p0, cfier );

    BOOST_REQUIRE_MESSAGE( rec.isEmpty(), "Expected that with base having only even element, and alt being empty, recombination should empty but is not." );
    BOOST_REQUIRE_MESSAGE( rec.isMatchAlt(), "Since alt is empty and recombination results in empty, rec should also match alt. However, this is not the case" );

    rec( c1_exp, c1_exp, cfier );
    BOOST_REQUIRE_MESSAGE( !rec.isEmpty(), "Expected that recombination of c1_exp with itself should not be empty." );
    BOOST_REQUIRE_MESSAGE( rec.isMatchBase(), "Expected that recombination of c1_exp with itself should match base");
    BOOST_REQUIRE_MESSAGE( rec.isMatchAlt(), "Expected that recombination  of c1_exp with itself should match alt");

    typename powerset_type::subset_ptr c = ps.create_subset( *rec.getResultSequence() );

    BOOST_REQUIRE_MESSAGE( *c == *c1_exp , "Expected that c == c1_exp");   
}

/**
 * Simple classifier test
 *
 * Create an empty powerset
 *
 * Create three empty subsets: P0, P1, C_exp;
 *
 * Alternate adding elements to P0 and P1.  Elements keys are normalized between [0, 1.0).
 *
 * C_exp is the subset with all elements.
 *
 * Recombine the sequences with P0 as the base sequence.  As a result of construction, P0 has
 * no odd index elements.  Therefore, any child sequence will inherit all 0s from P0.  Similarly,
 * P1 has no even elements.  Therefore, any child sequence will inherit all 0s from P1.  Hence,
 * this recombination results in the empty subset.
 *
 * Recombining the sequences with P1 as the base sequence results in the set with all elements
 * (C_exp).
 */
BOOST_AUTO_TEST_CASE( recombine_simple ) {
    powerset_type ps;

    typename powerset_type::subset_ptr p0 = ps.create_subset();
    typename powerset_type::subset_ptr p1 = ps.create_subset();
    typename powerset_type::subset_ptr c1_exp = ps.create_subset();

    const unsigned int sub_div = 3;
    const double div_offset = (1.0/ ((double)sub_div*bmap::bits_per_block));

    for( unsigned int i = 0; i < sub_div * bmap::bits_per_block; ++i ) {
        double v = (double) (i % bmap::bits_per_block);

        unsigned int s = (i / bmap::bits_per_block );
        
        double k = v / (double) bmap::bits_per_block + (double)s*div_offset;
        test_element te(k, 1.0, i);
        if( i & 1 ) {
            p1->addElement( te );
        } else {
            p0->addElement( te );
        }
        c1_exp->addElement(te);

        typename powerset_type::element_index_type idx = ps.find(te);

        BOOST_REQUIRE_MESSAGE( idx == i, "Unexpected index " << idx << " returned for " << i << "(" << k << ")" );
    }

    recombination_type rec;
    classifier_type cfier;

    rec( p0, p1, cfier );

    BOOST_REQUIRE_MESSAGE( rec.isEmpty(), "Expected that with p0 as base the simple recombination would be empty");

    rec( p1, p0, cfier );

    typename powerset_type::subset_ptr c = ps.create_subset( *rec.getResultSequence() );

    BOOST_REQUIRE_MESSAGE( *c == *c1_exp , "Expected that with p1 as base the simple recombination would be all set");   
}

/**
 * Range based recombination
 *
 * Create an empty powerset
 *
 * Create three empty subsets: P0, P1, C_exp;
 *
 * Define the upper bounds of recombination regions.
 *
 * Alternate adding elements to P0 and P1.  Elements keys are normalized between [0, 1.0).
 * Add elements to C_exp according to whether the parent subset it was added to, and
 * the value of its key.
 *
 * (-inf, 0.25) -> base
 * [0.25, 0.5) -> alt
 * [0.5, 0.75) -> base
 * [0.75, inf) -> alt
 *
 *
 */
BOOST_AUTO_TEST_CASE( recombine_range ) {
    powerset_type ps;

    typename powerset_type::subset_ptr p0 = ps.create_subset();
    typename powerset_type::subset_ptr p1 = ps.create_subset();
    typename powerset_type::subset_ptr c1_exp = ps.create_subset();
    typename powerset_type::subset_ptr c0_exp = ps.create_subset();

    // after recombination
    // [0, 0.25) -> p0; [0.25, 0.5) -> p1; [0.5, 0.75) -> p0; [0.75, 1.0) -> p1
    classifier2_type::region_upper_bounds pts(3);

    pts.push_back( test_element(0.25, 0.));
    pts.push_back( test_element(0.5, 0.));
    pts.push_back( test_element(0.75, 0.));

    const unsigned int sub_div = 3;
    const double div_offset = (1.0/ ((double)sub_div*bmap::bits_per_block));

    for( unsigned int i = 0; i < sub_div * bmap::bits_per_block; ++i ) {
        double v = (double) (i % bmap::bits_per_block);

        unsigned int s = (i / bmap::bits_per_block );
        
        double k = v / (double) bmap::bits_per_block + (double)s*div_offset;
        test_element te(k, 1.0);
        if( i & 1 ) {
            p1->addElement( te );
            if( (0.25 <= k && k < 0.5 ) ||
                (0.75 <= k && k < 1.0 ) ) {
                c1_exp->addElement(te);
            } else {
                c0_exp->addElement(te);
            }
        } else {
            p0->addElement( te );
            if( (0. <= k && k < 0.25 ) ||
                (0.5 <= k && k < 0.75) ) {
                c1_exp->addElement(te);
            } else {
                c0_exp->addElement(te);
            }
        }

        typename powerset_type::element_index_type idx = ps.find(te);

        BOOST_REQUIRE_MESSAGE( idx == i, "Unexpected index " << idx << " returned for " << i << "(" << k << ")" );
    }

    recombination2_type rec;
    classifier2_type cfier( pts );

    rec( p1, p0, cfier );

    typename powerset_type::subset_ptr c = ps.create_subset( *rec.getResultSequence() );

    BOOST_REQUIRE_MESSAGE( *c == *c1_exp, "Unexpected result: p1 + p0" );

    classifier2_type cfier2( pts );
    rec( p0, p1, cfier2 );
    typename powerset_type::subset_ptr c2 = ps.create_subset( *rec.getResultSequence() );

//    std::cerr << "P0:    " << *p0 << std::endl;
//    std::cerr << "P1:    " << *p1 << std::endl;
//    std::cerr << "C:     " << *c << std::endl;
//    std::cerr << "E(C1): " << *c1_exp << std::endl;
//    std::cerr << "C2:    " << *c2 << std::endl;
//    std::cerr << "E(C2): " << *c0_exp << std::endl;

    BOOST_REQUIRE_MESSAGE( *c2 == *c0_exp, "Unexpected result: p0 + p1" );
}
BOOST_AUTO_TEST_SUITE_END()
