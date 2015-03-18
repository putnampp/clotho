#ifndef SEQUENCE_MUTATOR_VARIABLE_HPP_
#define SEQUENCE_MUTATOR_VARIABLE_HPP_

#include "clotho/genetics/sequence_mutator_def.hpp"
#include "clotho/powerset/variable_subset.hpp"

template < class E, class B, class BM, class EK, class Generator >
class sequence_mutator< clotho::powersets::variable_subset< E, B, BM, EK >, Generator > {
public:
    typedef E           allele_type;
    typedef Generator   generator_type;

    typedef typename generator_type::result_type generator_result_type;

    typedef clotho::powersets::variable_subset< E, B, BM, EK > sequence_type;

    sequence_mutator( const generator_type & gen, unsigned int nEvents, unsigned int age = 0 ) : m_gen(gen), m_events(nEvents), m_age(age) {}

    void operator()( std::shared_ptr< sequence_type > & seq ) {
        operator()( *seq );
    }

    void operator()( sequence_type & seq ) {
        for( unsigned int i = 0; i < m_events; ++i ) {
            generator_result_type a = m_gen( seq, m_age );
            seq.addElement( a );
        }
    }

    unsigned int event_count() const {
        return m_events;
    }

protected:
    generator_type  m_gen;
    unsigned int    m_events, m_age;
};

#endif  // SEQUENCE_MUTATOR_VARIABLE_HPP_
