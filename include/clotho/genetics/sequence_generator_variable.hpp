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
#ifndef SEQUENCE_GENERATOR_VARIABLE_HPP_
#define SEQUENCE_GENERATOR_VARIABLE_HPP_

#include "clotho/genetics/sequence_generator_def.hpp"
#include "clotho/powerset/variable_subset.hpp"

template < class E, class B, class BM, class EK >
class sequence_generator< std::shared_ptr< clotho::powersets::variable_subset< E, B, BM, EK > > > {
public:
    typedef clotho::powersets::variable_subset< E, B, BM, EK >  sequence_type;
    typedef std::shared_ptr< sequence_type >                    result_type;
    typedef typename sequence_type::powerset_type               set_type;
    typedef sequence_generator< result_type >                   self_type;

    sequence_generator( set_type & s ) : m_set( &s ), m_res( s.create_subset() ) {}
    sequence_generator( const self_type & s ) : m_set( s.m_set ), m_res( s.m_res ) {}

    result_type operator()( ) {
//        return m_set->create_subset();
        return m_res;
    }

    virtual ~sequence_generator() {}
protected:
    set_type    * m_set;
    result_type m_res;
};

// sequence_generator< Sequence, sequence_generator< Sequence, recombine > >;
#endif  // SEQUENCE_GENERATOR_VARIABLE_HPP_
