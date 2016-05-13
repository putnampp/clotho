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
#ifndef CLOTHO_PHENOTYPE_EVALUATOR_HPP_
#define CLOTHO_PHENOTYPE_EVALUATOR_HPP_

#include <vector>

#include "clotho/utility/state_object.hpp"
#include "clotho/utility/log_helper.hpp"

namespace clotho {
namespace genetics {

template < class EvalMethodType >
class phenotype_evaluator {
public:
    typedef EvalMethodType                                              method_type;
    typedef typename accumulator_helper_of< method_type >::type         trait_accumulator_type;
    typedef typename accumulator_helper_of< method_type >::result_type  phenotype_type;

    typedef typename trait_accumulator_type::genetic_space_type     genetic_space_type;

    typedef typename trait_accumulator_type::trait_vector_type      trait_vector_type;
    typedef std::vector< phenotype_type >                           population_phenotype_type;

    typedef typename population_phenotype_type::iterator            phenotype_iterator;
    typedef typename population_phenotype_type::const_iterator      const_phenotype_iterator;

    typedef typename genetic_space_type::individual_iterator        individual_iterator;

    phenotype_evaluator( ) { }

    void update( genetic_space_type * gs, trait_accumulator_type & traits ) {
        size_t N = traits.size();
        size_t M = gs->individual_count();

        assert( M == N / 2 );

        resize( M );

        typedef typename genetic_space_type::individual_genome_type    individual_type;

        size_t i = 0;
        while( i < M )  {
            individual_type ind = gs->getIndividualAt( i );
            phenotype_type p = m_eval( traits.getTraitAt( ind.first ), traits.getTraitAt( ind.second ) );

            m_phenos[ i++ ] = p;
        }
    }

    phenotype_type & getPhenotypeAt( size_t idx ) {
        assert( 0 <= idx && idx < m_phenos.size() );

        return m_phenos[ idx ];
    }

    population_phenotype_type & getPhenotypes() {
        return m_phenos;
    }

    phenotype_iterator begin() {
        return m_phenos.begin();
    }

    const_phenotype_iterator begin() const {
        return m_phenos.begin();
    }

    phenotype_iterator end() {
        return m_phenos.end();
    }

    const_phenotype_iterator end() const {
        return m_phenos.end();
    }

    virtual ~phenotype_evaluator() {}

protected:

    void resize( size_t N ) {
        m_phenos.reserve( N );

        while( m_phenos.size() < N ) {
            m_phenos.push_back( phenotype_type() );
        }
    }

    population_phenotype_type   m_phenos;

    method_type      m_eval;
};

}   // namespace genetics
}   // namespace clotho

namespace clotho {
namespace utility {

template < class EvalMethodType >
struct state_getter< clotho::genetics::phenotype_evaluator< EvalMethodType > >  {
    typedef clotho::genetics::phenotype_evaluator< EvalMethodType > object_type;

    void operator()( boost::property_tree::ptree & s, object_type & obj ) {
        typedef typename object_type::phenotype_iterator iterator;

        iterator first = obj.begin(), last = obj.end();

        while( first != last ) {
            boost::property_tree::ptree p;
            clotho::utility::add_value_array( p, first->begin(), first->end() );
            ++first;

            s.push_back( std::make_pair( "", p ) );
        }
    }
};


}   // namespace utility
}   // namespace clotho

#endif  // CLOTHO_PHENOTYPE_EVALUATOR_HPP_
