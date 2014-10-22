#ifndef RECOMBINATION_DEF_HPP_
#define RECOMBINATION_DEF_HPP_

#include <cassert>

#include "clotho/powerset/variable_subset.hpp"

namespace clotho {
namespace recombine {

template < class Sequence >
class recombination {
public:
    typedef Sequence sequence_type;
    typedef Method method_type;

    void operator()( sequence_type base, sequence_type alt );
};

template < class Element, class Block, class BlockMap, class ElementKeyer >
class recombination< variable_subset< Element, Block, BlockMap, ElementKeyer > > {
public:
    typedef variable_subset< Element, Block, BlockMap, ElementKeyer >   subset_type;
    typedef typename subset_type::pointer                               sequence_type;

    typedef typename subset_type::bitset_type                           bit_sequence_type;
    typedef Block                                                       block_type;

    void operator()( sequence_type base, sequence_type alt ) {
        assert( base->isSameFamily( alt ) );
        m_match_base = m_match_alt = m_empty = true;

        if( base == alt ) {
            return;
        }

        m_res_seq.reset();
        m_res_seq.resize( base->max_size(), false );

        subset_type::cblock_iterator base_it = base->begin(), base_end = base->end();
        subset_type::cblock_iterator alt_it = base->begin(), alt_end = base->end();
        subset_type::block_iterator res_it = m_res_seq.m_bits.begin();

        while( true ) {
            if( alt_it == alt_end ) {
                while( base_it != base_end ) {
                    block_type _base = (*base_it++);
                    recombine( (*res_it), _base, 0 );
                    ++res_it;
                }
                break;
            }

            if( base_it == base_end ) {
                while( alt_it != alt_end ) {
                    block_type _alt = (*alt_it++);
                    recombine( (*res_it), 0, _alt );
                    ++res_it;
                }
                break;
            }

            block_type _base = (*base_it++), _alt = (*alt_it++);
            recombine( (*res_it), _base, _alt );
            ++res_it;
        }
    }

    bit_sequence_type * getResultSequence() { return &m_res_seq; }

    bool isMatchBase() const { return m_match_base; }
    bool isMatchAlt() const { return m_match_alt; }
    bool isEmpty() const    { return m_empty; }

    void recombine( block_type & res, const block_type base, const block_type alt ) {
        block_type hom = base & alt;
        block_type hets = base ^ alt;

        block_type bmask = base_mask_from_hets( hets );
        block_type amask = ((~bmask & hets) & alt);

        res =(hom | (amask | (bmask & base)));

        m_match_base = (m_match_base && (base == res) );
        m_match_alt = (m_match_alt && (alt == res) );
        m_empty = (m_empty && (res == 0) );
    }

    block_type base_mask_from_hets( block_type _bits ) {
        if( !_bits ) return (block_type)0;

        typename lowest_bit_map::block_type sub_block = (typename lowest_bit_map::block_type)(_bits);
        block_type res = 0;

        if( sub_block ) {    // block 0
            bit_walk( res, sub_block, 0 );
        }
        _bits /= lowest_bit_map::max_values;

        if( !_bits ) return res;

        sub_block = (typename lowest_bit_map::block_type)(_bits);

        if( sub_block ) {    // block 1
            bit_walk( res, sub_block, lowest_bit_map::block_width );
        }
        _bits /= lowest_bit_map::max_values;

        if( !_bits ) return res;

        sub_block = (typename lowest_bit_map::block_type)(_bits);

        if( sub_block ) {    // block 2
            bit_walk( res, sub_block, 2 * lowest_bit_map::block_width );
        }
        _bits /= lowest_bit_map::max_values;

        sub_block = (typename lowest_bit_map::block_type)(_bits);

        if( sub_block ) {    // block 3
            bit_walk( res, sub_block, 3 * lowest_bit_map::block_width );
        }

        return res;
    }

    inline void bit_walk( Block & res, typename lowest_bit_map::block_type sub_block, const Block & base, const Block & alt, set_iterator alpha_it, unsigned int _offset ) {
        const lowest_bit_map::value_type * v = low_bit_map.begin() + sub_block;
        do {
            unsigned int idx = _offset + v->bit_index;
            typename Alphabet::locus_t loc = accessor::get< set_iterator, typename Alphabet::locus_t >(alpha_it + idx);

            recombination_iterator rit = std::upper_bound( rec_points->begin(), rec_points->end(), loc);

            res |= ((((rit - rec_points->begin()) % 2) ? base : alt) & SortedAlleleAlphabet2::bit_position_masks[idx]);

            _offset += v->bit_shift_next;
            v = v->next_ptr;
        } while( v != NULL );
    }

    virtual ~recombination() {}
protected:
    bit_sequence_type   m_res_seq;
    bool    m_match_base, m_match_alt, m_empty;
};

}   // namespace recombine
}   // namespace clotho

#endif  // RECOMBINATION_DEF_HPP_
