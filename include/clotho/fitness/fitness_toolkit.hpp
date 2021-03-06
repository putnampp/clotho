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
#ifndef FITNESS_TOOLKIT_HPP_
#define FITNESS_TOOLKIT_HPP_

#include "clotho/fitness/ifitness_generator.hpp"

#include <unordered_map>
#include <memory>
#include <boost/property_tree/ptree.hpp>

#include <ostream>

class fitness_toolkit {
public:

    typedef std::unordered_map< std::string, ifitness_generator * >  generator_map;
    typedef generator_map::iterator                                  generator_iterator;
    typedef generator_map::const_iterator                            generator_citerator;

    static fitness_toolkit * getInstance() {
        static fitness_toolkit inst;
        return &inst;
    }

    std::shared_ptr< ifitness_generator > get_tool( boost::property_tree::ptree & config );

    void tool_configurations( boost::property_tree::ptree & config );

    void register_tool( ifitness_generator * gen );

    friend std::ostream & operator<<( std::ostream & out, const fitness_toolkit & ftk );

    virtual ~fitness_toolkit();

protected:
    fitness_toolkit();

    generator_map   m_tools;
};

std::ostream & operator<<( std::ostream & out, const fitness_toolkit & ftk );

#endif  // FITNESS_TOOLKIT_HPP_
