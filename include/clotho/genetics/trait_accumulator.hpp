#ifndef TRAIT_ACCUMULATOR_HPP_
#define TRAIT_ACCUMULATOR_HPP_

#include "trait_helper.hpp"
#include "clotho/utility/iterator_helper.hpp"
#include <vector>

template < class Element, class Weight >
class trait_accumulator {
public:
    typedef trait_accumulator< Element, Weight > self_type;
    typedef Element                         element_type;
    typedef typename Weight::value_type     value_type;
    typedef typename Weight::vector_type    result_type;

    typedef clotho::utility::iterator_helper< Element > element_helper;
    typedef typename element_helper::type               iterator;

    trait_accumulator( result_type & r) : m_res( &r ) {}
    trait_accumulator( const self_type & other ) : m_res( other.m_res ) {}

    void operator()( Element & elem ) {
        iterator first = element_helper::make_first(elem), last = element_helper::make_last(elem);

        typename result_type::iterator it = m_res->begin();
        while( first != last ) {
            if( it == m_res->end() ) {
                while( first != last ) {
                    m_res->push_back( *first );
                    ++first;
                }
                break;
            }
            (*it) += (*first);
            ++first;
            ++it;
        }
    }

    virtual ~trait_accumulator() {}
protected:
    result_type * m_res;
};

#endif  // TRAIT_ACCUMULATOR_HPP_
