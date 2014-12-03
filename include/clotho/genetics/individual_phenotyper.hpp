#ifndef INDIVIDUAL_PHENOTYPER_HPP_
#define INDIVIDUAL_PHENOTYPER_HPP_

#include "iterator_helper.hpp"

template < class TraitMatrix, class Environment >
class individual_phenotyper {
public:
    typedef individual_phenotyper< Genotyper, Environment > self_type;
    typedef Genotyper                                       genotyper_type;
    typedef Environment                                     environment_type;
    typedef Value                                           value_type;

    class result_type {
    public:
        friend class individual_phenotyper< Genotyper, Environment, Value >;

        
    };

    individual_phenotyper( genotyper_type & geno, environment_type & env ) :
        m_geno( geno )
        , m_env( env ) 
    {}

    template < class Individual >
    result_type operator()( Individual & ind ) {
        typedef iterator_helper< Individual > _iterator_helper;
        typedef typename _iterator_helper::type sequence_iterator;
        sequence_iterator first = iterator_helper::make_first(ind), last = iterator_helper::make_last( ind );

        result_type res;
        genotyper_type geno( res );

        std::for_each( first, last, geno );

        res += m_env(  );
        return res;
    }
    
protected:
    genotyper_type m_geno;
    environment_type m_env;
};

#endif  // INDIVIDUAL_PHENOTYPER_HPP_
