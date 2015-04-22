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
#ifndef POWERSET_HPP_
#define POWERSET_HPP_

#include "clotho/powerset/config.hpp"

#include <iostream>
#include <vector>
#include <unordered_map>
#include <set>
#include <boost/dynamic_bitset.hpp>

#include <cassert>
#include <limits>

#include <memory>

#include "clotho/utility/block_masks.hpp"
#include "clotho/powerset/element_key_of.hpp"
#include "clotho/powerset/block_map.hpp"

namespace clotho {
namespace powersets {

template < class Element, class Subset, class Block = unsigned long, class BlockMap = block_map< Element, Block > , class ElementKeyer = element_key_of< Element > >
class powerset {
public:
    typedef Element                 value_type;
    typedef Block                   block_type;

    typedef ElementKeyer                            element_keyer_type;
    typedef typename element_keyer_type::key_type   key_type;

    typedef std::vector< value_type >   set_type;
    typedef unsigned int                element_index_type;

    typedef BlockMap                block_map_type;

    typedef std::unordered_map< key_type, element_index_type > lookup_table_type;
    typedef typename lookup_table_type::iterator lookup_iterator;

    typedef typename set_type::iterator      fixed_iterator;
    typedef typename set_type::const_iterator cfixed_iterator;

    typedef typename set_type::iterator      variable_iterator;
    typedef typename set_type::const_iterator cvariable_iterator;

    typedef clotho::utility::block_masks< block_type > masks;

    static const unsigned int       max_elements = (std::numeric_limits< element_index_type >::max() / 2);
    //static const unsigned int       npos = -1;
    static const size_t             npos = -1;
    static const unsigned int       fixed_flag = ((element_index_type)npos ^ (((element_index_type)npos)>> 1));
    static const unsigned int       bits_per_block = masks::bits_per_block;

    typedef Subset                              subset_type;
    typedef std::shared_ptr< subset_type >            subset_ptr;
    typedef boost::dynamic_bitset< block_type > bitset_type;

    typedef std::set< subset_ptr >       family_type;
    typedef typename family_type::iterator  family_iterator;
    typedef typename family_type::const_iterator  cfamily_iterator;

    typedef block_type                           bitset_block_type;
    typedef typename bitset_type::allocator_type bitset_allocator_type;
    typedef std::vector< bitset_block_type, bitset_allocator_type > bitset_buffer_type;

    typedef typename bitset_buffer_type::iterator free_block_iterator;
    typedef typename bitset_buffer_type::const_iterator cfree_block_iterator;

    typedef typename subset_type::data_container_type dc_type;

    powerset();

    subset_ptr   empty_set() const;

    subset_ptr  create_subset();
    subset_ptr  create_subset( const dc_type & b );

//    subset_ptr  clone_subset( subset_ptr t );

//    void   copy_subset( subset_ptr t );
//    void   release_subset( subset_ptr t );

    /**
     *  Returns the element_index of the value, and if it was newly
     *
     *  bool - true if newly created, false if already exists
     */
    std::pair< element_index_type, bool > find_or_create( const value_type & v );

    element_index_type find_by_key( const key_type & k );
    element_index_type find( const value_type & v );

/// Element Iterators
    fixed_iterator          fixed_begin();
    cfixed_iterator    fixed_begin() const;

    fixed_iterator          fixed_end();
    cfixed_iterator    fixed_end() const;

    variable_iterator          variable_begin();
    cvariable_iterator    variable_begin() const;

    variable_iterator          variable_end();
    cvariable_iterator    variable_end() const;

    family_iterator     family_begin();
    cfamily_iterator    family_begin() const;

    family_iterator     family_end();
    cfamily_iterator    family_end() const;

    free_block_iterator free_begin();
    cfree_block_iterator free_begin() const;

    free_block_iterator free_end();
    cfree_block_iterator free_end() const;

    /**
     * Returns whether the set is empty
     */
    bool  empty() const;

    /**
     * Returns the number of keys in the lookup table
     */
    size_t size() const;

    /**
     *  Returns the number of elements which are variable in the set
     */
    size_t  variable_size() const;

    /**
     *  Returns the number of elements which are variable in the set
     */
    size_t  fixed_size() const;

    /**
     * Returns the number of variable elements which have been allocated
     */
    size_t  variable_allocated_size() const;

    /**
     *  Returns the number of free variable positions
     */
    size_t free_size() const;

    /**
     * Returns the number of subsets of the powerset exist
     */
    size_t family_size() const;

    bool isFamilyMember( subset_ptr p ) const;

    void clear();

/// Compact space
    void pruneSpace();

    virtual ~powerset();

//protected:
    CLOTHO_PROTECTED

    element_index_type addElement( const value_type & v );
    element_index_type appendElement( const value_type & v );

    element_index_type findFreeIndex( const value_type & v );

    unsigned int updateFixed();

    void updateFreeElements( subset_ptr s );
    void updateFreeElements( bitset_type & seq_data );
    void updateFreeElements( std::vector< element_index_type > & seq_data );

    void updateSubsetFree( subset_ptr s );
    void updateSubsetFree( std::vector< element_index_type > & seq_data );
    void updateSubsetFree( bitset_type & seq_data );

    void updateFreeIndex( element_index_type idx, bool state);

    void buildFreeRanges();

//    void clearGarbage();

    /**
     * Encode set index as a 1's complement
     *
     * Two element sets exist (fixed, variable).  An element may only exist in one set at any given time.
     * Variable elements have an index [0x0, 0x7FFFFFFF];
     * Fixed elements  have an index [0x80000000, 0xFFFFFFFE];
     * Index value of 0xFFFFFFFF is left available to indicate npos;
     *
     */
    inline element_index_type  encode_index( element_index_type idx, bool is_fixed ) {
        return ((is_fixed) ? (idx | fixed_flag) : idx);
    }

    /**
     * Decode set index
     *
     * Strips off high order bit of index
     *
     * The input index is modified if the value is fixed (negative).
     *
     * \return high order bit
     */
    inline bool decode_index( element_index_type & idx ) {
        if(idx & fixed_flag) {
            idx ^= fixed_flag;
            return true;
        }
        return false;
    }

    set_type            m_fixed, m_variable;
    lookup_table_type   m_lookup;


    family_type m_family;

//    unsigned int m_family_size; // running total number of active subset references

    /// META-INFORMATION about variable elements in powerset
    bitset_type     m_free_list, m_fixed_variable, m_lost_variable;
    bitset_type     m_variable_mask, m_free_ranges;

    block_map_type      m_block_map;
};

//
// Helper functions
//
template < class SubsetIterator >
void release_subsets( SubsetIterator first, SubsetIterator last ) {
    while( first != last ) {
        (*first)->release();
        ++first;
    }
}

//
// powerset Implementation Details
//

#define TEMPLATE_HEADER template < class Element, class Subset, class Block, class BlockMap, class ElementKeyer >
#define POWERSET_SPECIALIZATION powerset< Element, Subset, Block, BlockMap, ElementKeyer >

TEMPLATE_HEADER
POWERSET_SPECIALIZATION::powerset() //:
//    m_empty_set( this )
//    , m_family_size( 0 )
{}

TEMPLATE_HEADER
POWERSET_SPECIALIZATION::~powerset() {
    while( !m_family.empty() ) {
        subset_ptr s = *m_family.begin();
        m_family.erase( m_family.begin() );
        s.reset();
    }
}

TEMPLATE_HEADER
typename POWERSET_SPECIALIZATION::subset_ptr POWERSET_SPECIALIZATION::empty_set() const {
    return subset_ptr();
}

TEMPLATE_HEADER
typename POWERSET_SPECIALIZATION::subset_ptr POWERSET_SPECIALIZATION::create_subset( ) {
    subset_ptr sub( new subset_type( this ) );
    m_family.insert( sub );
    return sub;
}

TEMPLATE_HEADER
typename POWERSET_SPECIALIZATION::subset_ptr POWERSET_SPECIALIZATION::create_subset( const dc_type & b ) {
    subset_ptr sub( new subset_type( this, b ) );
    m_family.insert( sub );
    return sub;
}

//TEMPLATE_HEADER
//typename POWERSET_SPECIALIZATION::subset_ptr POWERSET_SPECIALIZATION::clone_subset( subset_ptr s ) {
//    subset_ptr sub( new subset_type( *s ));
//
//    m_family.insert( sub );
//    ++m_family_size;
//
//    return sub;
//}

//TEMPLATE_HEADER
//void POWERSET_SPECIALIZATION::copy_subset( subset_type * s ) {
//    assert( s->m_parent == this );
//    ++m_family_size;
//    ++s->m_ref_count;
//}

//TEMPLATE_HEADER
//void POWERSET_SPECIALIZATION::release_subset( subset_type * s ) {
//    assert( s->m_parent == this );
//    --m_family_size;
//    if( --s->m_ref_count == 0 ) {
//        m_family.erase( s );
//        m_garbage.push_back( s );
//    }
//}

TEMPLATE_HEADER
std::pair< typename POWERSET_SPECIALIZATION::element_index_type, bool > POWERSET_SPECIALIZATION::find_or_create( const value_type & v ) {
    std::pair< element_index_type, bool > res = std::make_pair( (element_index_type) npos, false );
    lookup_iterator it = m_lookup.find( element_keyer_type::get_key( v ) );

    if( it == m_lookup.end() ) {
        // DNE
        //
        res.first = addElement( v );
        res.second = true;
    } else {
        res.first = it->second;
    }
    return res;
}

TEMPLATE_HEADER
typename POWERSET_SPECIALIZATION::element_index_type POWERSET_SPECIALIZATION::find_by_key( const key_type & k ) {
    lookup_iterator it = m_lookup.find( k );

    return (( it == m_lookup.end() ) ? (element_index_type)npos : it->second);
}

TEMPLATE_HEADER
typename POWERSET_SPECIALIZATION::element_index_type POWERSET_SPECIALIZATION::find( const value_type & v ) {
    lookup_iterator it = m_lookup.find( element_keyer_type::get_key( v ) );

    return (( it == m_lookup.end() ) ? (element_index_type)npos : it->second);
}

TEMPLATE_HEADER
typename POWERSET_SPECIALIZATION::element_index_type POWERSET_SPECIALIZATION::addElement( const value_type & v ) {
    element_index_type idx = findFreeIndex(v);

    if( idx == (element_index_type) npos ) {
        idx = appendElement(v);
        m_lookup.insert( std::make_pair( element_keyer_type::get_key( v ), encode_index(idx, false)));
    } else {
        assert( idx < m_variable.size());

        if( m_variable[idx] != v ) {
            lookup_iterator it = m_lookup.find( element_keyer_type::get_key( m_variable[idx] ) );
            if( it != m_lookup.end() && it->second == idx ) {
                // if the element was lost
                m_lookup.erase( it );
            }

            m_lookup.insert( std::make_pair( element_keyer_type::get_key(v), encode_index(idx, false)));
        }
        m_variable[ idx ] = v;
    }

//    m_variable[ idx ] = v;
    assert( idx != (element_index_type) npos );

    updateFreeIndex( idx, false );

    return idx;
}

TEMPLATE_HEADER
typename POWERSET_SPECIALIZATION::element_index_type POWERSET_SPECIALIZATION::findFreeIndex( const value_type & v ) {
    element_index_type _offset = m_block_map( v );
    if( _offset == (element_index_type) npos ) {
        element_index_type idx = m_free_list.find_first();
        return idx;
    }

    block_type  mask = masks::position_mask( _offset );

    if( !m_free_ranges.empty() && (m_free_ranges.m_bits[0] & mask ) ) {
        if( m_free_list.num_blocks() == 1 ) {
            return _offset;
        }

        size_t path = 0, level = 2;
        while( level < m_free_ranges.num_blocks() ) {
            size_t p = (2 * path) + 1;

            if( !(m_free_ranges.m_bits[ p ] & mask ) ) {
                ++p;
                ++level;
            }
            path = p;
            level <<= 1;
        }

        size_t range_idx = (level - m_free_ranges.num_blocks());
        if(!( m_free_list.m_bits[ range_idx ] & mask) ) {
            ++range_idx;
        }
        _offset += range_idx * bits_per_block;
    } else {
        _offset = (element_index_type) npos;
    }

    return _offset;
}

TEMPLATE_HEADER
typename POWERSET_SPECIALIZATION::element_index_type POWERSET_SPECIALIZATION::appendElement( const value_type & v ) {
    element_index_type idx = m_variable.size();
    element_index_type _offset = m_block_map( v );

    idx += ((_offset != (element_index_type) npos ) ? _offset : 0 );

    unsigned int i = bits_per_block;
    while( i-- ) {
        m_variable.push_back( v );
    }

    m_free_list.resize( m_variable.size(), true );

    buildFreeRanges();

    return idx;
}

TEMPLATE_HEADER
void POWERSET_SPECIALIZATION::updateFreeIndex( element_index_type idx, bool state ) {
    m_free_list[idx] = state;

    element_index_type block_idx = idx / bits_per_block;//, block_offset = idx % bits_per_block;
    element_index_type path = block_idx + m_free_ranges.num_blocks();

//    block_type clear_mask = masks::position_mask( block_offset );
//   block_type set_mask = ((state) ? clear_mask : 0);

//    clear_mask = ~clear_mask;

//    m_free_list[idx] = state;
//    m_free_list.m_bits[ block_idx ] &= clear_mask;
//    m_free_list.m_bits[ block_idx ] |= set_mask;

    block_type c0 = m_free_list.m_bits[ block_idx ];

    block_idx += ((block_idx & 1) ? -1 : 1);

    block_type c1 = (( block_idx < m_free_list.num_blocks() ) ? m_free_list.m_bits[ block_idx ] : 0);

    if( path == 1 ) {
        m_free_ranges.m_bits[ 0 ] = (c0 | c1);
        return;
    }

    while( true ) {
        element_index_type p = (( path & 1 ) ? (path - 2) : (path - 1));
        p /= 2;
        path /= 2;

        m_free_ranges.m_bits[ p ] = (c0 | c1);

        if( !p ) break;

        c0 = m_free_ranges.m_bits[p];
        p += ((p & 1) ? 1 : -1);
        c1 = m_free_ranges.m_bits[p];
    }
}

TEMPLATE_HEADER
bool POWERSET_SPECIALIZATION::empty() const {
    return m_lookup.size() == 0;
}

TEMPLATE_HEADER
size_t POWERSET_SPECIALIZATION::size() const {
    return m_lookup.size();
}

TEMPLATE_HEADER
size_t POWERSET_SPECIALIZATION::variable_size() const {
    return m_variable.size() - m_free_list.count();
}

TEMPLATE_HEADER
size_t POWERSET_SPECIALIZATION::fixed_size() const {
    return m_fixed.size();
}

TEMPLATE_HEADER
size_t POWERSET_SPECIALIZATION::variable_allocated_size() const {
    return m_variable.size();
}

TEMPLATE_HEADER
size_t POWERSET_SPECIALIZATION::family_size() const {
    return m_family.size();
}

TEMPLATE_HEADER
bool POWERSET_SPECIALIZATION::isFamilyMember( subset_ptr p) const {
    return m_family.find(p) != m_family.end();
}

TEMPLATE_HEADER
size_t POWERSET_SPECIALIZATION::free_size() const {
    return m_free_list.count();
}

TEMPLATE_HEADER
void POWERSET_SPECIALIZATION::pruneSpace() {

    if( m_family.empty() ) return;

    if( m_fixed_variable.size() < m_variable.size() ) {
        m_fixed_variable.resize( m_variable.size(), false );
    }

    if( m_lost_variable.size() < m_variable.size() ) {
        m_lost_variable.resize( m_variable.size(), false );
    }

    m_fixed_variable.reset().flip();    // assume all elements are fixed
    m_lost_variable.reset();            // assume no elements are variable

    std::vector< family_iterator > to_remove;

    family_iterator it = m_family.begin();
    while( it != m_family.end() ) {
        if( (*it).use_count() > 1 ) {
            updateFreeElements( *it );
        } else {
            to_remove.push_back( it );
        }
        ++it;
    }

    m_lost_variable.flip();

    m_free_list = (m_fixed_variable | m_lost_variable);
    m_variable_mask = m_free_list;
    m_variable_mask.flip();         // not fixed and not lost

    // remove un-referenced subsets
    while( !to_remove.empty() ) {
        family_iterator tmp = to_remove.back();
        to_remove.pop_back();
        m_family.erase( tmp );
    }


    unsigned int fcount = updateFixed();
    if( fcount ) {
        it = m_family.begin();
        while( it != m_family.end() ) {
            updateSubsetFree( *it );
            ++it;
        }
    }

    buildFreeRanges();
}

TEMPLATE_HEADER
void POWERSET_SPECIALIZATION::updateSubsetFree( subset_ptr s ) {
    updateSubsetFree( s->m_data );
}

TEMPLATE_HEADER
void POWERSET_SPECIALIZATION::updateSubsetFree( bitset_type & seq_data ) {
    typedef typename bitset_buffer_type::iterator       buffer_iterator;
    typedef typename bitset_buffer_type::const_iterator buffer_citerator;

    assert( seq_data.num_blocks() <= m_variable_mask.num_blocks() );

    buffer_citerator fit = m_variable_mask.m_bits.begin();
    buffer_iterator  sit = seq_data.m_bits.begin(), send = seq_data.m_bits.end();

    while( sit != send ) {
        (*sit) &= (*fit);
        ++sit;
        ++fit;
    }
}

TEMPLATE_HEADER
void POWERSET_SPECIALIZATION::updateSubsetFree( std::vector< element_index_type > & seq_data ) {

    typedef std::vector< element_index_type >::iterator data_iterator;

    data_iterator res = seq_data.begin(), first = seq_data.begin(), last = seq_data.end();

    while( first != last ) {
        if( m_variable_mask[ *first ] ) {
            *res = *first;
            ++res;
        }
        ++first;
    }

    seq_data.erase( res, seq_data.end() );
}

TEMPLATE_HEADER
void POWERSET_SPECIALIZATION::updateFreeElements( subset_ptr s ) {
    updateFreeElements( s->m_data );
}

TEMPLATE_HEADER
void POWERSET_SPECIALIZATION::updateFreeElements( bitset_type & seq_data ) {
    typedef typename bitset_buffer_type::iterator       buffer_iterator;
    typedef typename bitset_buffer_type::const_iterator buffer_citerator;

    // don't need b/c empty_set not stored in family
    //if( s == empty_set() ) return;

    buffer_citerator cit = seq_data.m_bits.begin(), cend = seq_data.m_bits.end();

    assert( seq_data.num_blocks() <= m_fixed_variable.num_blocks() );

    buffer_iterator int_it = m_fixed_variable.m_bits.begin(), int_end = m_fixed_variable.m_bits.end();
    buffer_iterator un_it = m_lost_variable.m_bits.begin();

    // handle common prefix
    while( cit != cend ) {
        block_type b = (*cit);
        ++cit;

        (*int_it) &= b;
        ++int_it;

        (*un_it) |= b;
        ++un_it;
    }

    // handle suffix
    while( int_it != int_end ) {
        (*int_it) = 0;
        ++int_it;
    }
}

TEMPLATE_HEADER
void POWERSET_SPECIALIZATION::updateFreeElements( std::vector< element_index_type > & seq_data ) {
    if( seq_data.empty() ) {
        m_fixed_variable.reset();
        return;
    }

    typedef std::vector< element_index_type >::const_iterator data_citerator;
    data_citerator it = seq_data.begin();

    bitset_type _set( m_variable.size(), 0 );

    // explode the index set into a bit map
    do {
        _set[ *it ] = true;
    } while( ++it != seq_data.end() );

    m_fixed_variable &= _set;
    m_lost_variable |= _set;
}

TEMPLATE_HEADER
typename POWERSET_SPECIALIZATION::fixed_iterator POWERSET_SPECIALIZATION::fixed_begin() {
    return m_fixed.begin();
}

TEMPLATE_HEADER
typename POWERSET_SPECIALIZATION::cfixed_iterator POWERSET_SPECIALIZATION::fixed_begin() const {
    return m_fixed.begin();
}

TEMPLATE_HEADER
typename POWERSET_SPECIALIZATION::fixed_iterator POWERSET_SPECIALIZATION::fixed_end() {
    return m_fixed.end();
}

TEMPLATE_HEADER
typename POWERSET_SPECIALIZATION::cfixed_iterator POWERSET_SPECIALIZATION::fixed_end() const {
    return m_fixed.end();
}

TEMPLATE_HEADER
typename POWERSET_SPECIALIZATION::variable_iterator POWERSET_SPECIALIZATION::variable_begin() {
    return m_variable.begin();
}

TEMPLATE_HEADER
typename POWERSET_SPECIALIZATION::cvariable_iterator POWERSET_SPECIALIZATION::variable_begin() const {
    return m_variable.begin();
}

TEMPLATE_HEADER
typename POWERSET_SPECIALIZATION::variable_iterator POWERSET_SPECIALIZATION::variable_end() {
    return m_variable.end();
}

TEMPLATE_HEADER
typename POWERSET_SPECIALIZATION::cvariable_iterator POWERSET_SPECIALIZATION::variable_end() const {
    return m_variable.end();
}

TEMPLATE_HEADER
typename POWERSET_SPECIALIZATION::family_iterator POWERSET_SPECIALIZATION::family_begin() {
    return m_family.begin();
}

TEMPLATE_HEADER
typename POWERSET_SPECIALIZATION::cfamily_iterator POWERSET_SPECIALIZATION::family_begin() const {
    return m_family.begin();
}

TEMPLATE_HEADER
typename POWERSET_SPECIALIZATION::family_iterator POWERSET_SPECIALIZATION::family_end() {
    return m_family.end();
}

TEMPLATE_HEADER
typename POWERSET_SPECIALIZATION::cfamily_iterator POWERSET_SPECIALIZATION::family_end() const {
    return m_family.end();
}

TEMPLATE_HEADER
typename POWERSET_SPECIALIZATION::free_block_iterator POWERSET_SPECIALIZATION::free_begin() {
    return m_free_list.m_bits.begin();
}

TEMPLATE_HEADER
typename POWERSET_SPECIALIZATION::cfree_block_iterator POWERSET_SPECIALIZATION::free_begin() const {
    return m_free_list.m_bits.begin();
}

TEMPLATE_HEADER
typename POWERSET_SPECIALIZATION::free_block_iterator POWERSET_SPECIALIZATION::free_end() {
    return m_free_list.m_bits.end();
}

TEMPLATE_HEADER
typename POWERSET_SPECIALIZATION::cfree_block_iterator POWERSET_SPECIALIZATION::free_end() const {
    return m_free_list.m_bits.end();
}

TEMPLATE_HEADER
unsigned int POWERSET_SPECIALIZATION::updateFixed() {
    typename bitset_type::size_type idx = m_fixed_variable.find_first();

    unsigned int fcount = 0;
    while( idx != bitset_type::npos ) {
        element_index_type eid = m_fixed.size();
        m_fixed.push_back( m_variable[ idx ] );

        lookup_iterator it = m_lookup.find( element_keyer_type::get_key( m_fixed.back() ) );
        assert( it != m_lookup.end() );

        it->second = encode_index( eid, true );
        idx = m_fixed_variable.find_next(idx);
        ++fcount;
    }

    return fcount;
}

TEMPLATE_HEADER
void POWERSET_SPECIALIZATION::buildFreeRanges() {
    m_free_ranges.reset();

    if( m_variable.size() == 0 || m_free_list.num_blocks() == 0 ) return;

    if( m_free_list.size() < m_variable.size() ) {
        m_free_list.resize( m_variable.size(), false );
    }

    size_t inode_count = 1;
    size_t bcount = (m_variable.size() / bits_per_block ) + 1;

    while( inode_count < bcount ) {
        inode_count <<= 1;
    }

    m_free_ranges.resize( inode_count * bits_per_block, false );

    assert( m_free_ranges.num_blocks() >= m_free_list.num_blocks() );

    if( m_free_list.num_blocks() <= 2 ) {
        m_free_ranges.m_bits[0] = m_free_list.m_bits[0];
        if( m_free_list.num_blocks() == 2 )
            m_free_ranges.m_bits[0] |= m_free_list.m_bits[1];

        return;
    }

    size_t b_idx = m_free_list.num_blocks() - 1,
           r_idx = ((b_idx + m_free_ranges.num_blocks()) / 2 ) - 1;

    if( !(b_idx & 1 ) ) {
        // even index; (odd number of blocks)
        m_free_ranges.m_bits[ r_idx-- ] = m_free_list.m_bits[b_idx--];
    }

    while( b_idx < m_free_list.num_blocks() ) {
        m_free_ranges.m_bits[ r_idx ] = m_free_list.m_bits[ b_idx-- ];
        m_free_ranges.m_bits[ r_idx--] |= m_free_list.m_bits[ b_idx-- ];
    }

    b_idx = m_free_ranges.num_blocks() - 2;

    while( r_idx < m_free_ranges.num_blocks() ) {
        m_free_ranges.m_bits[ r_idx ] = m_free_ranges.m_bits[ b_idx-- ];
        m_free_ranges.m_bits[ r_idx-- ] |= m_free_ranges.m_bits[ b_idx-- ];
    }
}

}   // namespace powersets
}   // namespace clotho

#undef POWERSET_SPECIALIZATION
#undef TEMPLATE_HEADER

#ifdef UNDEF_FRIEND
#undef BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS
#undef UNDEF_FRIEND
#endif  // UNDEF_FRIEND

#endif  // POWERSET_HPP_
