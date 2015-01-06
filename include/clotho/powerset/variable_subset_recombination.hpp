#ifndef CLOTHO_VARIABLE_SUBSET_RECOMBINATION_DEF_HPP_
#define CLOTHO_VARIABLE_SUBSET_RECOMBINATION_DEF_HPP_

#include "clotho/powerset/config.hpp"

#include <iostream>
#include <cassert>
#include <algorithm>

#include "clotho/recombination/recombination_def.hpp"
#include "clotho/recombination/bit_block_recombiner.hpp"

#include "clotho/powerset/variable_subset.hpp"

namespace clotho {
namespace recombine {

template < class Element, class Block, class BlockMap, class ElementKeyer, class Classifier >
class recombination< clotho::powersets::variable_subset< Element, Block, BlockMap, ElementKeyer >, Classifier > {
public:
    typedef clotho::powersets::variable_subset< Element, Block, BlockMap, ElementKeyer >   subset_type;
    typedef typename subset_type::pointer                               sequence_type;

    typedef Classifier                                                  classifier_type;

    typedef clotho::recombine::bit_block_recombiner< Classifier >        recombiner_type;

    typedef typename subset_type::bitset_type                           bit_sequence_type;
    typedef Block                                                       block_type;

    static const unsigned int bits_per_block = sizeof( block_type ) * 8;

    recombination( ) :
        m_classifier()
        , m_swap(false) {
    }

    recombination( const classifier_type & c, bool sw = false ) :
        m_classifier( c )
        , m_swap( sw ) {
    }

    void operator()( sequence_type base, sequence_type alt ) {
        if( m_swap ) {
            execute( alt, base, m_classifier );

            std::swap( m_match_base, m_match_alt );
        } else {
            execute( base, alt, m_classifier );
        }
    }

    void operator()( sequence_type base, sequence_type alt, classifier_type & elem_classifier ) {
        execute( base, alt, elem_classifier );
    }

    void execute( sequence_type base, sequence_type alt, classifier_type & elem_classifier ) {
        m_match_base = m_match_alt = m_empty = true;
        m_res_seq.reset();

        if( base == alt ) {
            // pointers match
            if( base ) {
                m_empty = false;
                m_res_seq.resize( base->max_size(), false);
                std::copy( base->begin(), base->end(), m_res_seq.m_bits.begin() );
            }
            return;
        }

        typename subset_type::cblock_iterator base_it, base_end;
        typename subset_type::cblock_iterator alt_it, alt_end;

        typename subset_type::powerset_type::cvariable_iterator elem_it, elem_end;

        if( !base ) {
            // base sequence is empty
            // NOTE: since base != alt, alt must be defined
            alt_it = alt->begin();
            alt_end = alt->end();

            // slight "trick" to cause following loop to iterate only
            // over the alt sequence
            base_it = alt->end();
            base_end = alt->end();

            elem_it = alt->getParent()->variable_begin();
            elem_end = alt->getParent()->variable_end();

            m_res_seq.resize( alt->max_size(), false );
        } else if( !alt ) {
            // alt sequence is empty
            base_it = base->begin();
            base_end = base->end();

            // slight "trick" to cause following loop to iterate only
            // over the base sequence
            alt_it = base->end();
            alt_end = base->end();

            elem_it = base->getParent()->variable_begin();
            elem_end = base->getParent()->variable_end();

            m_res_seq.resize( base->max_size(), false );
        } else {
            assert( base->isSameFamily( alt ) );

            base_it = base->begin();
            base_end = base->end();

            alt_it = alt->begin();
            alt_end = alt->end();

            elem_it = base->getParent()->variable_begin();
            elem_end = base->getParent()->variable_end();

            m_res_seq.resize( base->max_size(), false );
        }

        typename subset_type::block_iterator res_it = m_res_seq.m_bits.begin();

        recombiner_type  brecombiner( elem_classifier );

        while( true ) {
            if( alt_it == alt_end ) {
                while( base_it != base_end ) {
                    assert( elem_it != elem_end );
                    block_type _base = (*base_it++);
                    block_type r = brecombiner( _base, (block_type)0, elem_it );

                    m_match_base = (m_match_base && (_base == r) );
                    m_match_alt = (m_match_alt && ((block_type)0 == r) );
                    m_empty = (m_empty && (r == (block_type)0) );

                    (*res_it++) = r;
                    elem_it += bits_per_block;
                }
                break;
            }

            if( base_it == base_end ) {
                while( alt_it != alt_end ) {
                    assert( elem_it != elem_end );
                    block_type _alt = (*alt_it++);
                    block_type r = brecombiner( (block_type)0, _alt, elem_it );

                    m_match_base = (m_match_base && ((block_type)0 == r) );
                    m_match_alt = (m_match_alt && (_alt == r) );
                    m_empty = (m_empty && (r == (block_type)0) );

                    (*res_it++) = r;
                    elem_it += bits_per_block;
                }
                break;
            }

            assert( elem_it != elem_end );
            block_type _base = (*base_it++), _alt = (*alt_it++);
            block_type r = brecombiner( _base, _alt, elem_it );

            m_match_base = (m_match_base && (_base == r) );
            m_match_alt = (m_match_alt && (_alt == r) );
            m_empty = (m_empty && (r == 0) );

            (*res_it++) = r;
            elem_it += bits_per_block;
        }
    }

    bit_sequence_type * getResultSequence() {
        return &m_res_seq;
    }

    bool isMatchBase() const {
        return m_match_base;
    }
    bool isMatchAlt() const {
        return m_match_alt;
    }
    bool isEmpty() const    {
        return m_empty;
    }

    virtual ~recombination() {}
protected:
    classifier_type     m_classifier;
    bool                m_swap;
    bit_sequence_type   m_res_seq;
    bool    m_match_base, m_match_alt, m_empty;
};

}   // namespace recombine
}   // namespace clotho

#endif  // CLOTHO_VARIABLE_SUBSET_RECOMBINATION_DEF_HPP_
