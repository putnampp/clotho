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
#ifndef LINEAR_POPULATION_GROWTH_HPP_
#define LINEAR_POPULATION_GROWTH_HPP_

#include "clotho/genetics/ipopulation_growth.hpp"

extern const std::string LINEAR_POP_NAME;

/**
 * Linear Population Growth
 *
 * Integer truncation of float value
 *
 * size = (unsigned int) (a * psize + b)
 */
class linear_population_growth : public ipopulation_growth {
public:
    linear_population_growth( double a = 1.0, double b = 0.0 );

    unsigned int operator()( unsigned int psize, unsigned int );

    const std::string name() const;

    void log( std::ostream & out ) const;

    virtual ~linear_population_growth();

protected:
    double m_A, m_B;
};

#endif  // LINEAR_POPULATION_GROWTH_HPP_
