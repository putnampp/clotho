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
#ifndef INDIVIDUAL_GENERATOR_HPP_
#define INDIVIDUAL_GENERATOR_HPP_

template < class Population, class SelectionModel, class ReproductionModel >
class individual_generator;

#include "individual_selector.hpp"
#include "individual_reproduction.hpp"

//template < class IndividualType, class URNG, class MutationModel, class RecombinationModel, class ReproMethodTag >
//class individual_generator< std::vector< IndividualType >, individual_selector< URNG >, individual_reproduction< IndividualType, MutationModel, RecombinationModel, ReproMethodTag > > {
template < class IndividualType, class Selector, class MutationModel, class RecombinationModel, class ReproMethodTag >
class individual_generator< std::vector< IndividualType >, Selector, individual_reproduction< IndividualType, MutationModel, RecombinationModel, ReproMethodTag > > {
public:
    typedef std::vector< IndividualType >   population_type;
    typedef IndividualType                  individual_type;
    typedef IndividualType                  result_type;

//    typedef individual_selector< URNG >     selection_type;
    typedef Selector                        selection_type;
    typedef individual_reproduction< IndividualType, MutationModel, RecombinationModel, ReproMethodTag > reproduction_type;

    individual_generator( population_type * pop, selection_type & sel, reproduction_type & repro, unsigned int age = 0 ) :
        m_pop( pop )
        , m_sel( sel )
        , m_repro( repro )
        , m_gen( age ) {
    }

    result_type operator()() {
        individual_type p0 = m_pop->at(m_sel()), p1 = m_pop->at( m_sel());

        return m_repro(p0, p1, m_gen);
    }
protected:
    population_type * m_pop;
    selection_type  m_sel;
    reproduction_type m_repro;
    unsigned int    m_gen;
};


#endif  // INDIVIDUAL_GENERATOR_HPP_
