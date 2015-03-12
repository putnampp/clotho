#ifndef CLOTHO_VECTOR_SUBSET_RECOMBINATION_DEF_HPP_
#define CLOTHO_VECTOR_SUBSET_RECOMBINATION_DEF_HPP_

#include "clotho/powerset/config.hpp"

#include <iostream>
#include <cassert>
#include <algorithm>

#include "clotho/recombination/recombination_def.hpp"
#include "clotho/recombination/bit_block_recombiner.hpp"

#include "clotho/powerset/vector_subset.hpp"

#include "clotho/recombination/inspect_methods.hpp"

namespace clotho {
namespace recombine {

template < class Element, class Block, class BlockMap, class ElementKeyer, class Classifier, class T1 >
class recombination< clotho::powersets::vector_subset< Element, Block, BlockMap, ElementKeyer >, Classifier, clotho::recombine::inspection::tag::copy_matching_classify_mismatch, T1 > {
public:
    typedef clotho::powersets::vector_subset< Element, Block, BlockMap, ElementKeyer >   subset_type;
    typedef typename subset_type::pointer                               sequence_type;

    typedef Classifier                                                  classifier_type;

//    typedef typename subset_type::bitset_type                           bit_sequence_type;
    typedef typename subset_type::data_container_type                   data_container_type;
    typedef typename subset_type::element_index_type                    element_index_type;

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
        m_res_seq.clear();

        if( base == alt ) {
            // pointers match
            if( base ) {
                m_res_seq.reserve(base->size());
                std::copy( base->begin(), base->end(), m_res_seq.begin() );
                m_empty = m_res_seq.empty();
            }
            return;
        }

        typename subset_type::element_index_iterator base_it, base_end;
        typename subset_type::element_index_iterator alt_it, alt_end;

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
        } else {
            assert( base->isSameFamily( alt ) );

            base_it = base->begin();
            base_end = base->end();

            alt_it = alt->begin();
            alt_end = alt->end();

            elem_it = base->getParent()->variable_begin();
            elem_end = base->getParent()->variable_end();
        }

        while( true ) {
            if( alt_it == alt_end ) {
                while( base_it != base_end ) {
                    element_index_type _base = (*base_it);
                    ++base_it;

                    if( elem_classifier( *(elem_it + _base ) ) ) {
                        m_res_seq.push_back( _base );
                        m_match_alt = false;
                    } else {
                        m_match_base = false;
                    }
                }
                break;
            }

            if( base_it == base_end ) {
                while( alt_it != alt_end ) {
                    element_index_type _alt = (*alt_it);
                    ++alt_it;

                    if( !elem_classifier( *(elem_it + _alt) ) ) {
                        m_res_seq.push_back( _alt );
                        m_match_base = false;
                    } else {
                        m_match_alt = false;
                    }
                }
                break;
            }

            element_index_type _base = (*base_it), _alt = (*alt_it);

            if( _base < _alt ) {
                if( elem_classifier( *(elem_it + _base) ) ) {
                    m_res_seq.push_back(_base);
                    m_match_alt = false;
                } else {
                    m_match_base = false;
                }
                ++base_it;
            } else if( _alt < _base ) {
                if( ! elem_classifier( *(elem_it + _alt) ) ) {
                    m_res_seq.push_back(_alt);
                    m_match_base = false;
                } else {
                    m_match_alt = false;
                }
                ++alt_it;
            } else {
                m_res_seq.push_back(_base);
                ++base_it;
                ++alt_it;
            }
        }

        m_empty = m_res_seq.empty();
    }

    data_container_type * getResultSequence() {
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
    data_container_type m_res_seq;
    bool    m_match_base, m_match_alt, m_empty;
};

}   // namespace recombine
}   // namespace clotho

#endif  // CLOTHO_VECTOR_SUBSET_RECOMBINATION_DEF_HPP_