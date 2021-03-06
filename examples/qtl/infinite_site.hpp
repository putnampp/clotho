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
#ifndef INFINITE_SITE_GENERATOR_HPP_
#define INFINITE_SITE_GENERATOR_HPP_

template < class Sequence >
struct infinite_site {
    typedef Sequence type;
};

#include "clotho/powerset/variable_subset.hpp"
#include "clotho/powerset/vector_subset.hpp"
#include "clotho/powerset/powerset_no_dup_pred.hpp"
#include "clotho/utility/random_generator.hpp"

namespace clotho {
namespace utility {

template < class URNG, class E, class B, class BM, class EK >
class random_generator< URNG, infinite_site< clotho::powersets::variable_subset< E, B, BM, EK > > > {
public:
    typedef E result_type;
    typedef clotho::powersets::variable_subset< E, B, BM, EK > sequence_type;

    typedef clotho::mutations::no_duplicate_pred< typename sequence_type::powerset_type > no_dup_type;
    typedef clotho::utility::random_generator< URNG, E >    random_element;

    random_generator( URNG & rng, boost::property_tree::ptree & config ) :
        m_rng( &rng )
        , m_gen( rng, config ) {
    }

    result_type operator()( std::shared_ptr< sequence_type > & seq, unsigned int age = 0 ) {
        no_dup_type tester( seq->getParent() );
        return generate( tester, age );
    }

    result_type operator()( const sequence_type & seq, unsigned int age = 0 ) {
        no_dup_type tester( seq.getParent() );
        return generate( tester, age );
    }

protected:
    inline result_type generate( no_dup_type & tester, unsigned int age ) {
        result_type res = m_gen(age);
        while( !tester(res) ) {
            res = m_gen(age);
        }
        return res;
    }

    URNG * m_rng;
    random_element m_gen;
};

template < class URNG, class E, class B, class BM, class EK >
class random_generator< URNG, infinite_site< clotho::powersets::vector_subset< E, B, BM, EK > > > {
public:
    typedef E result_type;
    typedef clotho::powersets::vector_subset< E, B, BM, EK > sequence_type;

    typedef clotho::mutations::no_duplicate_pred< typename sequence_type::powerset_type > no_dup_type;
    typedef clotho::utility::random_generator< URNG, E >    random_element;

    random_generator( URNG & rng, boost::property_tree::ptree & config ) :
        m_rng( &rng )
        , m_gen( rng, config ) {
    }

    result_type operator()( std::shared_ptr< sequence_type > & seq, unsigned int age = 0 ) {
        no_dup_type tester(seq->getParent());
        return generate(tester, age);
    }

    result_type operator()( const sequence_type & seq, unsigned int age = 0 ) {
        no_dup_type tester(seq.getParent());
        return generate( tester, age);
    }

protected:

    inline result_type generate( no_dup_type & tester, unsigned int age ) {
        result_type res = m_gen(age);
        while( !tester(res) ) {
            res = m_gen(age);
        }
        return res;
    }
    URNG * m_rng;
    random_element m_gen;
};
}   // namespace utility
}   // namespace clotho

#endif  // INFINITE_SITE_HPP_
