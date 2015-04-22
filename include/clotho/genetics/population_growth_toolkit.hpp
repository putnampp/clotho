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
#ifndef POPULATION_GROWTH_TOOLKIT_HPP_
#define POPULATION_GROWTH_TOOLKIT_HPP_

#include "clotho/genetics/ipopulation_growth_generator.hpp"

#include <unordered_map>
#include <memory>
#include <boost/property_tree/ptree.hpp>

#include <ostream>

extern const std::string POPULATION_GROWTH_BLOCK_K;

class population_growth_toolkit {
public:

    typedef std::unordered_map< std::string, ipopulation_growth_generator * >  generator_map;
    typedef generator_map::iterator                                  generator_iterator;
    typedef generator_map::const_iterator                            generator_citerator;

    static population_growth_toolkit * getInstance() {
        static population_growth_toolkit inst;
        return &inst;
    }

    std::shared_ptr< ipopulation_growth_generator > get_tool( boost::property_tree::ptree & config );

    void tool_configurations( boost::property_tree::ptree & config );

    void register_tool( ipopulation_growth_generator * gen );

    friend std::ostream & operator<<( std::ostream & out, const population_growth_toolkit & ftk );

    virtual ~population_growth_toolkit();

protected:
    population_growth_toolkit();

    generator_map   m_tools;
};

std::ostream & operator<<( std::ostream & out, const population_growth_toolkit & ftk );

#endif  // POPULATION_GROWTH_TOOLKIT_HPP_
