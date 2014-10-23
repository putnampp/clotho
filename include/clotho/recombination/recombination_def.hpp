#ifndef RECOMBINATION_DEF_HPP_
#define RECOMBINATION_DEF_HPP_

#include <cassert>

#include "clotho/powerset/variable_subset.hpp"
#include "clotho/utility/bit_walker.hpp"

namespace clotho {
namespace recombine {

template < class Sequence, class Classifier >
class recombination;

template < class Element, class Block, class BlockMap, class ElementKeyer, class Classifier >
class recombination< clotho::powersets::variable_subset< Element, Block, BlockMap, ElementKeyer >, Classifier > {
public:
    typedef clotho::powersets::variable_subset< Element, Block, BlockMap, ElementKeyer >   subset_type;
    typedef typename subset_type::pointer                               sequence_type;

    typedef Classifier                                                  classifier_type;

    typedef typename subset_type::bitset_type                           bit_sequence_type;
    typedef Block                                                       block_type;

    typedef clotho::utility::block_walker< block_type, unsigned short > block_walker_type;

    void operator()( sequence_type base, sequence_type alt, classifier_type & elem_classifier) {
        assert( base->isSameFamily( alt ) );
        m_match_base = m_match_alt = m_empty = true;

        if( base == alt ) {
            return;
        }

        m_res_seq.reset();
        m_res_seq.resize( base->max_size(), false );

        typename subset_type::cblock_iterator base_it = base->begin(), base_end = base->end();
        typename subset_type::cblock_iterator alt_it = base->begin(), alt_end = base->end();
        typename subset_type::block_iterator res_it = m_res_seq.m_bits.begin();

        while( true ) {
            if( alt_it == alt_end ) {
                while( base_it != base_end ) {
                    block_type _base = (*base_it++);
                    recombine( (*res_it), _base, 0, elem_classifier );
                    ++res_it;
                }
                break;
            }

            if( base_it == base_end ) {
                while( alt_it != alt_end ) {
                    block_type _alt = (*alt_it++);
                    recombine( (*res_it), 0, _alt, elem_classifier );
                    ++res_it;
                }
                break;
            }

            block_type _base = (*base_it++), _alt = (*alt_it++);
            recombine( (*res_it), _base, _alt, elem_classifier );
            ++res_it;
        }
    }

    bit_sequence_type * getResultSequence() { return &m_res_seq; }

    bool isMatchBase() const { return m_match_base; }
    bool isMatchAlt() const { return m_match_alt; }
    bool isEmpty() const    { return m_empty; }

    void recombine( block_type & res, const block_type base, const block_type alt, classifier_type & elem_classifier ) {
        block_type hom = base & alt;
        block_type hets = base ^ alt;

        elem_classifier.resetResult();

        block_walker_type::apply( hets, elem_classifier );

        block_type base_mask = elem_classifier.getResult();
        block_type alt_mask = ((~base_mask & hets) & alt);

        res =(hom | (alt_mask | (base_mask & base)));

        m_match_base = (m_match_base && (base == res) );
        m_match_alt = (m_match_alt && (alt == res) );
        m_empty = (m_empty && (res == 0) );
    }

    virtual ~recombination() {}
protected:
    bit_sequence_type   m_res_seq;
    bool    m_match_base, m_match_alt, m_empty;
};

}   // namespace recombine
}   // namespace clotho

#endif  // RECOMBINATION_DEF_HPP_
