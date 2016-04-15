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
#ifndef RECOMBINER_HPP_
#define RECOMBINER_HPP_

template < class Sequence, class Engine >
class recombiner {
public:
    template < class Generator >
    Sequence operator()( Sequence base, Sequence alt, Generator & gen, bool should_copy ) {
        assert( false );
    }
};

#include "clotho/powerset/variable_subset_recombination.hpp"

template < class E, class B, class BM, class EK, class C, class T0, class T1 >
class recombiner< clotho::powersets::variable_subset< E, B, BM, EK >,
    clotho::recombine::recombination< clotho::powersets::variable_subset< E, B, BM, EK >, C, T0, T1 > > {
public:
    typedef clotho::powersets::variable_subset< E, B, BM, EK >  sequence_type;
    typedef typename sequence_type::pointer                     sequence_pointer;

    typedef sequence_pointer                                    result_type;

    typedef C                                                   classifier_type;
//    typedef typename C::param_type                              classifier_param;

    typedef clotho::recombine::recombination< clotho::powersets::variable_subset< E, B, BM, EK >, C, T0, T1 > engine_type;

    recombiner( const engine_type & eng ) : m_engine( eng ) {}

    result_type operator()( std::pair< sequence_pointer, sequence_pointer > ind, bool should_copy ) {
        return operator()( ind.first, ind.second, should_copy );
    }

    result_type operator()( sequence_pointer base, sequence_pointer alt, bool should_copy ) {
        if( base == alt ) {
            if( !base ) {
                // both sequences must be NULL
                return sequence_pointer();
            } else if( should_copy ) {
                return base;
            }

            sequence_pointer res = base->clone();
            return res;
        }

        // at this point neither base or alt can both be NULL
        typename sequence_type::powerset_type * p = ((base) ? base->getParent() : alt->getParent());
        assert( p );

        m_engine( base, alt );

        sequence_pointer res;
        if( !m_engine.isEmpty() ) {
            if( !should_copy ) {
                res = p->create_subset( *m_engine.getResultSequence() );
            } else if( m_engine.isMatchBase() ) {
                res = base;
            } else if( m_engine.isMatchAlt() ) {
                res = alt;
            } else {
                res = p->create_subset( *m_engine.getResultSequence() );
            }
        } else {
            res = p->create_subset();
        }

        return res;
    }

protected:
    engine_type   m_engine;
};

#include "clotho/powerset/vector_subset_recombination.hpp"

template < class E, class B, class BM, class EK, class C, class T0, class T1 >
class recombiner< clotho::powersets::vector_subset< E, B, BM, EK >,
    clotho::recombine::recombination< clotho::powersets::vector_subset< E, B, BM, EK >, C, T0, T1 > > {
public:
    typedef clotho::powersets::vector_subset< E, B, BM, EK >  sequence_type;
    typedef typename sequence_type::pointer                     sequence_pointer;

    typedef sequence_pointer                                    result_type;

    typedef C                                                   classifier_type;

    typedef clotho::recombine::recombination< clotho::powersets::vector_subset< E, B, BM, EK >, C, T0, T1 > engine_type;

    recombiner( const engine_type & eng ) : m_engine( eng ) {}

    result_type operator()( std::pair< sequence_pointer, sequence_pointer > ind, bool should_copy ) {
        return operator()( ind.first, ind.second, should_copy );
    }

    result_type operator()( sequence_pointer base, sequence_pointer alt, bool should_copy ) {
        if( base == alt ) {
            if( !base ) {
                // both sequences must be NULL
                return sequence_pointer();
            } else if( should_copy ) {
                return base;
            }

            sequence_pointer res = base->clone();
            return res;
        }

        // at this point neiter base or alt can both be NULL
        typename sequence_type::powerset_type * p = ((base) ? base->getParent() : alt->getParent());
        assert( p );

        m_engine( base, alt );

        sequence_pointer res;
        if( !m_engine.isEmpty() ) {
            if( !should_copy ) {
                res = p->create_subset( *m_engine.getResultSequence() );
            } else if( m_engine.isMatchBase() ) {
                res = base;
            } else if( m_engine.isMatchAlt() ) {
                res = alt;
            } else {
                res = p->create_subset( *m_engine.getResultSequence() );
            }
        } else {
            res = p->create_subset();
        }

        return res;
    }

protected:
    engine_type   m_engine;
};

#include "recombination_generator.hpp"
#include "recombiner_generator.hpp"

#endif  // RECOMBINER_HPP_
