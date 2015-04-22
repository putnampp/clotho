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

struct basic_allele {
    double  m_key;
    double  m_select, m_dom;

    bool    m_neutral;

    basic_allele( double k, double sel = 0.0, double dom = 1.0, bool neut = true ) :
        m_key(k), m_select(sel), m_dom(dom), m_neutral( neut ) {
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
    typedef double key_type;

    static key_type get_key( const basic_allele & e ) {
        return e.m_key;
    }
};

template <>
struct normalized_key< basic_allele > : public key_range < > {
    typedef double key_type;

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
void fitness_method< double, multiplicative_heterozygous_tag >::operator()< basic_allele >( double & res, const basic_allele & elem, double scale ) {
    res *= (1. + elem.m_select * elem.m_dom);
}

template < > template < >
void fitness_method< double, multiplicative_homozygous_tag >::operator()< basic_allele >( double & res, const basic_allele & elem, double scale ) {
    res *= (1. + elem.m_select * scale);
}

template < > template < >
void fitness_method< double, additive_heterozygous_tag >::operator()< basic_allele >( double & res, const basic_allele & elem, double scale ) {
    res += (elem.m_select * elem.m_dom);
}

template < > template < >
void fitness_method< double, additive_homozygous_tag >::operator()< basic_allele >( double & res, const basic_allele & elem, double scale ) {
    res += (elem.m_select * scale);
}

}   // namespace fitness {
}   // namespace clotho {

#endif  // BASIC_ALLELE_H_
