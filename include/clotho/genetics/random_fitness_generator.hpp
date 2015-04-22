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
#ifndef RANDOM_FITNESS_GENERATOR_HPP_
#define RANDOM_FITNESS_GENERATOR_HPP_

#include "clotho/genetics/ifitness_generator.hpp"
#include "clotho/genetics/random_fitness_metric.hpp"

/**
 * Computes the phenotype scaling factor based upon:
 * -  the current population's size (N)
 * -  population mutation rate (mu)
 * -  user-defined scaling factor (s)
 *
 * Computes the theoretical standard deviation of the
 * population, and scales it accordingly.
 *
 * Generates a random fitness metric.
 */
class random_fitness_generator : public ifitness_generator {
public:
    typedef random_fitness_metric result_type;

    random_fitness_generator();
    random_fitness_generator( boost::property_tree::ptree & config );

    std::shared_ptr< ifitness_generator > create( boost::property_tree::ptree & config ) const;

    std::shared_ptr< ifitness > generate( const std::vector< std::vector< double > > & pop_traits );

    const std::string name() const;

    void log( std::ostream & out ) const;

    virtual ~random_fitness_generator();

protected:

    void parseConfig( boost::property_tree::ptree & config );

    double m_scale, m_mu;
};

#endif  // RANDOM_FITNESS_GENERATOR_HPP_
