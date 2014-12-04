#ifndef INDIVIDUAL_PHENOTYPER_HPP_
#define INDIVIDUAL_PHENOTYPER_HPP_

//#include "iterator_helper.hpp"
//
#include "individual_trait_accumulator.hpp"

template < class IndividualType, class Environment >
class individual_phenotyper {
public:
    typedef individual_phenotyper< IndividualType, Environment > self_type;
    typedef IndividualType                                  individual_type;
    typedef individual_trait_accumulator< IndividualType >  individual_trait_acc_type;
    typedef Environment                                     environment_type;

    typedef typename individual_trait_acc_type::result_type result_type;

    individual_phenotyper( environment_type & env ) : m_env(env) {}

    result_type operator()( individual_type & ind ) {
        individual_trait_acc_type acc;
        result_type res = acc( ind );

        typename result_type::iterator it = res.begin();
        while( it != res.end() ) {
            (*it) += m_env( );
            ++it;
        }

        return res;
    }
protected:
    environment_type    m_env;
};

struct no_type {};

template < class IndividualType >
class individual_phenotyper< IndividualType, no_type > {
public:
    typedef individual_phenotyper< IndividualType, no_type > self_type;
    typedef IndividualType                                  individual_type;
    typedef individual_trait_accumulator< IndividualType >  individual_trait_acc_type;
    typedef no_type                                         environment_type;

    typedef typename individual_trait_acc_type::result_type result_type;

    individual_phenotyper( ) {}

    result_type operator()( individual_type & ind ) {
        individual_trait_acc_type acc;
        result_type res = acc( ind );

        return res;
    }
};

#endif  // INDIVIDUAL_PHENOTYPER_HPP_
