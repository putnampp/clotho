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
#ifndef INDIVIDUAL_PHENOTYPER_HPP_
#define INDIVIDUAL_PHENOTYPER_HPP_

#include <boost/property_tree/ptree.hpp>
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

    individual_phenotyper( boost::property_tree::ptree & config ) : m_env(config) {}
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
    individual_phenotyper( boost::property_tree::ptree & config ) {}

    result_type operator()( individual_type & ind ) {
        individual_trait_acc_type acc;
        result_type res = acc( ind );

        return res;
    }
};

#endif  // INDIVIDUAL_PHENOTYPER_HPP_
