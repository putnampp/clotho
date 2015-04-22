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
#ifndef PREDEFINED_POPULATION_GROWTH_HPP_
#define PREDEFINED_POPULATION_GROWTH_HPP_

#include "clotho/genetics/ipopulation_growth.hpp"

#include <vector>

extern const std::string PREDEF_POP_NAME;

/**
 * Predefined Population Growth
 *
 * Population Sizes are cyclically returned from a predefined list of
 * per generation sizes.
 */
class predefined_population_growth : public ipopulation_growth {
public:
    predefined_population_growth( const std::vector< unsigned int > & sizes );

/**
 * Population Size is returned for the generation
 *
 * The previous population size is ignored
 */
    unsigned int operator()( unsigned int , unsigned int generation );

    const std::string name() const;

    void log( std::ostream & out ) const;

    virtual ~predefined_population_growth();

protected:
    std::vector< unsigned int > m_gen_size;
};

#endif  // PREDEFINED_POPULATION_GROWTH_HPP_
