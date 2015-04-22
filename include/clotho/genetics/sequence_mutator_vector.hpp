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
#ifndef SEQUENCE_MUTATOR_VECTOR_HPP_
#define SEQUENCE_MUTATOR_VECTOR_HPP_

#include "clotho/genetics/sequence_mutator_def.hpp"
#include "clotho/powerset/vector_subset.hpp"

template < class E, class B, class BM, class EK, class Generator >
class sequence_mutator< clotho::powersets::vector_subset< E, B, BM, EK >, Generator > {
public:
    typedef E           allele_type;
    typedef Generator   generator_type;

    typedef typename generator_type::result_type generator_result_type;

    typedef clotho::powersets::vector_subset< E, B, BM, EK > sequence_type;

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
    unsigned int    m_events;
    unsigned int    m_age;
};

#endif  // SEQUENCE_MUTATOR_VECTOR_HPP_
