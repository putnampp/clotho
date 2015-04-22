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

#include "test_element.h"
#include "clotho/utility/bit_block_iterator.hpp"
#include "clotho/powerset/variable_subset.hpp"
#include "clotho/powerset/variable_subset_recombination.hpp"

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
template < class Block >
struct odd_from_base_classifier {
    typedef Block block_type;

    odd_from_base_classifier() : m_res(0) {}

/**
 * Does nothing
 */
    inline void updateOffset( unsigned int idx ) { return; }

/**
 * Classifies all indices which are odd as being (1 - base)
 */
    void operator()( unsigned int idx ) {
        if( idx & 1 ) {
            m_res |= ((block_type)1 << (idx % (sizeof(block_type) * 8)));
        }
    }

/**
 * Resets the result to be all 0s
 */
    void resetResult() {
        m_res = (block_type)0;
    }

/**
 * Returns the most recent classified indices 
 */
    block_type getResult() const {
        return m_res;
    }

    block_type m_res;
};

typedef odd_from_base_classifier< Block > classifier_type;
typedef clotho::recombine::recombination< subset_type, classifier_type > recombination_type;

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
template < class Block >
struct range_classifier {
    typedef Block block_type;
    typedef typename powerset_type::element_keyer_type keyer_type;
    typedef std::vector< keyer_type::key_type > recombination_points;
    typedef recombination_points::iterator      rec_iterator;

    range_classifier( const recombination_points & pts, typename powerset_type::cvariable_iterator vit) :
        m_res(0), m_pts(pts), m_vit(vit) {
        std::sort( m_pts.begin(), m_pts.end() );
    }

/**
 * Classifies a element at (m_vit + idx) as being either (1 - base) or not (0 - alt)
 * 
 * Sets corresponding bit at idx of result accordingly
 */
    void operator()( unsigned int idx ) {
        keyer_type::key_type k = m_keyer( *(m_vit + idx) );

        rec_iterator rit = std::upper_bound( m_pts.begin(), m_pts.end(), k );

        if( (rit - m_pts.begin()) % 2 == 0 ) {
            m_res |= ((block_type)1 << idx);
        }
    }

/**
 * Used to advance the m_vit iterator by offset elements
 */
    inline void updateOffset( unsigned int offset ) {
        m_vit += offset;
    }

/**
 * Resets the result classification to all 0
 */
    void resetResult() {
        m_res = (block_type)0;
    }

/**
 * Returns the result of the most recent classified indices since last reset
 */
    block_type getResult() const {
        return m_res;
    }

    block_type              m_res;
    recombination_points    m_pts;
    typename powerset_type::cvariable_iterator m_vit;
    keyer_type              m_keyer;
};

typedef range_classifier< Block > classifier2_type;
typedef clotho::recombine::recombination< subset_type, classifier2_type > recombination2_type;

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
    //

    recombination_type rec;
    classifier_type cfier;

    rec( p0, p1, cfier );

    BOOST_REQUIRE_MESSAGE( rec.isEmpty(), "Expected that with p0 and p1 being empty, recombination would be empty");
    BOOST_REQUIRE_MESSAGE( rec.isMatchBase(), "Expected that with p0 and p1 being empty, recombination should match base");
    BOOST_REQUIRE_MESSAGE( rec.isMatchAlt(), "Expected that with p0 and p1 being empty, recombination should match alt");

    test_element te( 0.5, 1.0 );
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
        test_element te(k, 1.0);
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
    classifier2_type::recombination_points pts(3);

    pts.push_back(0.25);
    pts.push_back(0.5);
    pts.push_back(0.75);

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
    classifier2_type cfier( pts, ps.variable_begin() );

    rec( p1, p0, cfier );

    typename powerset_type::subset_ptr c = ps.create_subset( *rec.getResultSequence() );

    BOOST_REQUIRE_MESSAGE( *c == *c1_exp, "Unexpected result: p1 + p0" );

    classifier2_type cfier2( pts, ps.variable_begin() );
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
