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
#ifndef BASIC_ALLELE_H_
#define BASIC_ALLELE_H_

#include <ostream>

static const double DEFAULT_SELECTION = 0.0;
static const double DEFAULT_DOMINANCE = 1.0;

#ifndef ALL_NEUTRAL_OPTIMIZATION
static const bool   DEFAULT_NEUTRAL = true;
#else
static const bool   DEFAULT_NEUTRAL = false;
#endif

extern const string ALLELE_BLOCK_K;

struct basic_allele {
    typedef double real_type;
    typedef double key_type;

    key_type    m_key;
    real_type   m_select, m_dom;

    bool    m_neutral;
    unsigned int m_age;

    basic_allele( key_type k= 0.0, real_type sel = 0.0, real_type dom = 1.0, bool neut = true, unsigned int age = 0 ) :
        m_key(k), m_select(sel), m_dom(dom), m_neutral( neut ), m_age(age) {
    }

    bool isNeutral() const {
        return m_neutral;
    }

    virtual ~basic_allele() {}
};

inline bool operator<( const basic_allele & lhs, const basic_allele & rhs ) {
    return lhs.m_key < rhs.m_key;
}

inline bool operator==( const basic_allele & lhs, const basic_allele & rhs ) {
    return (lhs.m_key == rhs.m_key);
}

inline bool operator!=( const basic_allele & lhs, const basic_allele & rhs ) {
    return (lhs.m_key != rhs.m_key);
}

std::ostream & operator<<( std::ostream & lhs, const basic_allele & rhs ) {
    lhs << "{" << rhs.m_key << ";" << rhs.m_select << ";" << rhs.m_dom << "}";
    return lhs;
}

#include "clotho/powerset/element_key_of.hpp"
#include "clotho/powerset/normalized_key.hpp"

namespace clotho {
namespace powersets {

template <>
struct element_key_of< basic_allele > {
    typedef typename basic_allele::key_type key_type;

    static key_type get_key( const basic_allele & e ) {
        return e.m_key;
    }
};

template <>
struct normalized_key< basic_allele > : public key_range < > {
    typedef typename basic_allele::key_type key_type;

    static key_type get_key( const basic_allele & a ) {
        return a.m_key;
    }
};

}   // namespace powersets {
}   // namespace clotho {

#include "clotho/fitness/fitness_method.hpp"

namespace clotho {
namespace fitness {

template < > template < >
void fitness_method< typename basic_allele::real_type, multiplicative_heterozygous_tag >::operator()< basic_allele >( typename basic_allele::real_type & res, const basic_allele & elem, typename basic_allele::real_type scale ) {
    res *= (1. + elem.m_select * elem.m_dom);
}

template < > template < >
void fitness_method< typename basic_allele::real_type, multiplicative_homozygous_tag >::operator()< basic_allele >( typename basic_allele::real_type & res, const basic_allele & elem, typename basic_allele::real_type scale ) {
    res *= (1. + elem.m_select * scale);
}

template < > template < >
void fitness_method< typename basic_allele::real_type, additive_heterozygous_tag >::operator()< basic_allele >( typename basic_allele::real_type & res, const basic_allele & elem, typename basic_allele::real_type scale ) {
    res += (elem.m_select * elem.m_dom);
}

template < > template < >
void fitness_method< typename basic_allele::real_type, additive_homozygous_tag >::operator()< basic_allele >( typename basic_allele::real_type & res, const basic_allele & elem, typename basic_allele::real_type scale ) {
    res += (elem.m_select * scale);
}

}   // namespace fitness {
}   // namespace clotho {

#include "basic_allele_generator.hpp"

#include "neutral_method.hpp"

template <>
bool neutral_method::test< basic_allele >( const basic_allele & ba ) {
    return ba.isNeutral();
}

#endif  // BASIC_ALLELE_H_
