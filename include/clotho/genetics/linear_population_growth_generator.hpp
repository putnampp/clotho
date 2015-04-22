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
#ifndef LINEAR_POPULATION_GROWTH_GENERATOR_HPP_
#define LINEAR_POPULATION_GROWTH_GENERATOR_HPP_

#include "clotho/genetics/ipopulation_growth_generator.hpp"
#include "clotho/genetics/linear_population_growth.hpp"

class linear_population_growth_generator : public ipopulation_growth_generator {
public:
    typedef linear_population_growth result_type;

    linear_population_growth_generator();
    linear_population_growth_generator( boost::property_tree::ptree & config );

    std::shared_ptr< ipopulation_growth_generator > create( boost::property_tree::ptree & config ) const;

    std::shared_ptr< ipopulation_growth > generate() const;

    const std::string name() const;

    void log( std::ostream & out ) const;

    virtual ~linear_population_growth_generator();

protected:
    void parseConfig( boost::property_tree::ptree & config );

    double m_A, m_B;
};

#endif  // LINEAR_POPULATION_GROWTH_GENERATOR_HPP_
