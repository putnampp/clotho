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
#include "clotho/fitness/constant_fitness_generator.hpp"
#include "clotho/fitness/fitness_toolkit.hpp"
#include "clotho/utility/clotho_strings.hpp"

constant_fitness_generator::constant_fitness_generator() {
    fitness_toolkit::getInstance()->register_tool( this );
}

constant_fitness_generator::constant_fitness_generator( boost::property_tree::ptree & config ) {}

std::shared_ptr< ifitness_generator > constant_fitness_generator::create( boost::property_tree::ptree & config ) const {
    std::shared_ptr< ifitness_generator > t( new constant_fitness_generator( config ) );
    return t;
}

std::shared_ptr< ifitness > constant_fitness_generator::generate( const std::vector< std::vector< real_type > > & pop_traits ) {
    std::shared_ptr< ifitness > t( new constant_fitness_metric( 1.0 ) );
    return t;
}

const std::string constant_fitness_generator::name() const {
    return CONSTANT_NAME;
}

void constant_fitness_generator::log( std::ostream & out ) const {
    out << "{" << CONSTANT_NAME << "}\n";
}

constant_fitness_generator::~constant_fitness_generator() {}

static const constant_fitness_generator cfg_reg;
