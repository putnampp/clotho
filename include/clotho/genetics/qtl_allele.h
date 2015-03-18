#ifndef QTL_ALLELE_H_
#define QTL_ALLELE_H_

#include "basic_allele.h"
#include "trait_weight.hpp"

struct qtl_allele : public basic_allele {

    typedef trait_weight< double >              weight_type;
    typedef typename weight_type::vector_type   trait_weights;
    typedef typename trait_weights::iterator    weight_iterator;
    typedef typename trait_weights::const_iterator weight_citerator;

    qtl_allele( double k, double sel, double dom, bool neut, unsigned int age, const trait_weights & coeff ) :
        basic_allele( k, sel, dom, neut, age ), m_weights( coeff ) {
    }

    qtl_allele( const basic_allele & all, const trait_weights & coeff ) :
        basic_allele( all ), m_weights( coeff ) {
    }

    weight_iterator    begin() {
        return m_weights.begin();
    }
    weight_citerator   begin() const {
        return m_weights.begin();
    }

    weight_iterator    end() {
        return m_weights.end();
    }
    weight_citerator   end() const {
        return m_weights.end();
    }

    virtual ~qtl_allele() {}

    trait_weights   m_weights;
};

inline bool operator<( const qtl_allele & lhs, const qtl_allele & rhs ) {
    return lhs.m_key < rhs.m_key;
}

inline bool operator==( const qtl_allele & lhs, const qtl_allele & rhs ) {
    return lhs.m_key == rhs.m_key;
}

inline bool operator!=( const qtl_allele & lhs, const qtl_allele & rhs ) {
    return lhs.m_key != rhs.m_key;
}

inline std::ostream & operator<<( std::ostream & os, const qtl_allele & qtl ) {
    os << "{" << qtl.m_key << ";" << qtl.m_select << ";" << qtl.m_dom << ";<";
    qtl_allele::weight_citerator it = qtl.begin();
    if( it != qtl.end() ) {
        os << *it;
        while( ++it != qtl.end() ) {
            os << "," << *it;
        }
    }
    os << ">;" << qtl.m_age << "}";
    return os;
}

#include "clotho/powerset/element_key_of.hpp"
#include "clotho/powerset/normalized_key.hpp"

namespace clotho {
namespace powersets {

template <>
struct element_key_of< qtl_allele > {
    typedef double key_type;

    static key_type get_key( const qtl_allele & e ) {
        return e.m_key;
    }
};

template <>
struct normalized_key< qtl_allele > : public key_range < > {
    typedef double key_type;

    static key_type get_key( const qtl_allele & a ) {
        return a.m_key;
    }
};

}   // namespace powersets {
}   // namespace clotho {

#include "clotho/fitness/fitness_method.hpp"

namespace clotho {
namespace fitness {

template < > template < >
void fitness_method< double, multiplicative_heterozygous_tag >::operator()< qtl_allele >( double & res, const qtl_allele & elem, double scale ) {
    res *= (1. + elem.m_select * elem.m_dom);
}

template < > template < >
void fitness_method< double, multiplicative_homozygous_tag >::operator()< qtl_allele >( double & res, const qtl_allele & elem, double scale ) {
    res *= (1. + elem.m_select * scale);
}

template < > template < >
void fitness_method< double, additive_heterozygous_tag >::operator()< qtl_allele >( double & res, const qtl_allele & elem, double scale ) {
    res += (elem.m_select * elem.m_dom);
}

template < > template < >
void fitness_method< double, additive_homozygous_tag >::operator()< qtl_allele >( double & res, const qtl_allele & elem, double scale ) {
    res += (elem.m_select * scale);
}

}   // namespace fitness {
}   // namespace clotho {

#include "qtl_allele_generator.hpp"
#include "qtl_region_classifier.hpp"
#include "qtl_weight_iterator.hpp"

#include "neutral_method.hpp"

template < >
bool neutral_method::test< qtl_allele >( const qtl_allele & all ) {
    return all.isNeutral();
}

#endif  // QTL_ALLELE_H_
