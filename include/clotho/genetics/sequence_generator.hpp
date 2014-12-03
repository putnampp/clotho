#ifndef SEQUENCE_GENERATOR_HPP_
#define SEQUENCE_GENERATOR_HPP_

template < class Sequence >
class sequence_generator;

#include "clotho/powerset/variable_subset.hpp"

template < class E, class B, class BM, class EK >
class sequence_generator< std::shared_ptr< clotho::powersets::variable_subset< E, B, BM, EK > > > {
public:
    typedef clotho::powersets::variable_subset< E, B, BM, EK >  sequence_type;
    typedef std::shared_ptr< sequence_type >                    result_type;
    typedef typename sequence_type::powerset_type               set_type;
    typedef sequence_generator< result_type >                   self_type;

    sequence_generator( set_type & s ) : m_set( &s ) {}
    sequence_generator( const self_type & s ) : m_set( s.m_set ) {}

    result_type operator()( ) {
        return m_set->create_subset();
    }

    virtual ~sequence_generator() {}
protected:
    set_type    * m_set;
};

// sequence_generator< Sequence, sequence_generator< Sequence, recombine > >;
#endif  // SEQUENCE_GENERATOR_HPP_
