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
#ifndef POPULATION_MUTATION_GENERATOR_HPP_
#define POPULATION_MUTATION_GENERATOR_HPP_

#include <boost/property_tree/ptree.hpp>
#include "clotho/cuda/mutation/mutation_event_generator.hpp"

template < class PopulationType >
class population_mutation_generator {
public:
    typedef PopulationType                              population_type;

    typedef typename population_type::allele_space_type allele_space_type;
    typedef MutationEventGenerator< allele_space_type >        generator_type;

    population_mutation_generator( boost::property_tree::ptree & config ) :
        m_gen(config)
    {}

    void operator()( allele_space_type * space, unsigned int N) {

    }

    virtual ~population_mutation_generator() {}

protected:
    generator_type  m_gen;
};

#endif  // POPULATION_MUTATION_GENERATOR_HPP_
