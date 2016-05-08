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
#ifndef CONSTANT_FITNESS_GENERATOR_HPP_
#define CONSTANT_FITNESS_GENERATOR_HPP_

#include "clotho/fitness/ifitness_generator.hpp"
#include "clotho/fitness/constant_fitness_metric.hpp"

class constant_fitness_generator : public ifitness_generator {
public:
    typedef double                  real_type;
    typedef constant_fitness_metric result_type;

    constant_fitness_generator();
    constant_fitness_generator( boost::property_tree::ptree & config );

    std::shared_ptr< ifitness_generator > create( boost::property_tree::ptree & config ) const;

    std::shared_ptr< ifitness > generate( const std::vector< std::vector< real_type > > & pop_traits );
    std::shared_ptr< ifitness > generate( size_t );

    const std::string name() const;

    void log( std::ostream & out ) const;

    virtual ~constant_fitness_generator();

protected:

    void parseConfig( boost::property_tree::ptree & config );

    double m_scale;
};

#endif  // CONSTANT_FITNESS_GENERATOR_HPP_
