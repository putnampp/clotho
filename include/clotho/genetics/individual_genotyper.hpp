#ifndef INDIVIDUAL_GENOTYPER_HPP_
#define INDIVIDUAL_GENOTYPER_HPP_

#include <vector>
#include <algorithm>

#include "iterator_helper.hpp"

struct add {
    template < class ResultType >
    ResultType operator()( ResultType a, ResultType b ) { return a + b; }
};

template < class SetType, class ValueType = double, class OP = add >
class individual_genotyper {
public:
    typedef individual_genotyper< SetType, ResultType > self_type;
    typedef ValueType       value_type;
    typedef SetType         set_type;

    typedef typename set_type::iterator iterator;

    class result_type {
    public:
        typedef std::vector< value_type > values_type;
        typedef typename values_type::iterator iterator;
        typedef typename values_type::const_iterator citerator;

        friend class individual_genotyper< SetType, ValueType, OP >;

        result_type ( size_t n = 1 ) : m_values( n ) {}
        result_type ( const result_type & r ) : m_values( r.m_values ) {}

        iterator    begin() { return m_values.begin(); }
        citerator   begin() const { return m_values.begin(); }

        iterator    end() { return m_values.end(); }
        citerator   end() const { return m_values.end(); }

    protected:
        values_type m_values;
    };

    individual_genotyper( set_type & s, size_t n = 1 ) : m_set( &s ), m_nTraits(n) {}

    template < class Sequence >
    result_type operator()( Sequence & seq ) {

        typedef iterator_helper< Sequence >      allele_helper;
        typedef typename allele_helper::type     allele_iterator;

        allele_iterator first = allele_helper::make_first( seq ), last = allele_helper::make_last( seq );

        std::for_each( first, last, weight );
    }

protected:
    set_type    * m_set;
    size_t      m_nTraits;
};

#endif  // INDIVIDUAL_GENOTYPER_HPP_
