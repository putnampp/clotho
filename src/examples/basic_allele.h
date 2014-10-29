#ifndef BASIC_ALLELE_H_
#define BASIC_ALLELE_H_

#include <ostream>

struct basic_allele {
    double m_key;
    double m_select, m_dom;

    basic_allele( double k, double sel = 0.0, double dom = 1.0 ) :
        m_key(k), m_select(sel), m_dom(dom)
    {}

    bool isNeutral() const { return (m_select == 0.0); }

    virtual ~basic_allele() {}
};

bool operator==( const basic_allele & lhs, const basic_allele & rhs ) {
    return (lhs.m_key == rhs.m_key);
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

    static key_type get_key( const basic_allele & e ) { return e.m_key; }
};

template <>
struct normalized_key< basic_allele > : public key_range < > {
    typedef double key_type;

    static key_type get_key( const basic_allele & a ) { return a.m_key; }
};

}   // namespace powersets {
}   // namespace clotho {

#include "clotho/mutation/element_generator.hpp"

#endif  // BASIC_ALLELE_H_