#ifndef INDIVIDUAL_INITIALIZER_HPP_
#define INDIVIDUAL_INITIALIZER_HPP_

template < class Individual, class Generator >
class individual_initializer;

#include "clotho/powerset/variable_subset.hpp"
#include "sequence_generator.hpp"

template < class Sequence >
class individual_initializer< std::pair< Sequence, Sequence >, sequence_generator< Sequence > > {
public:
    typedef sequence_generator< Sequence >      generator_type;
    typedef std::pair< Sequence, Sequence >     result_type;

    individual_initializer( const generator_type & gen ) : m_gen( gen ) {}

    result_type operator()() {
        return std::make_pair( m_gen(), m_gen() );
    }
protected:
    generator_type m_gen;
};

#endif  // INDIVIDUAL_INITIALIZER_HPP_
