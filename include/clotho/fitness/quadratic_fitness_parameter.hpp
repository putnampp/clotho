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
#ifndef QUADRATIC_FITNESS_PARAMETER_HPP_
#define QUADRATIC_FITNESS_PARAMETER_HPP_

#include <boost/property_tree/ptree.hpp>
#include "clotho/utility/clotho_strings.hpp"

extern const std::string QUAD_NAME;

template < class RealType >
struct quadratic_fitness_parameter {
    RealType    m_scale;

    quadratic_fitness_parameter( RealType s = 1.0) : m_scale( s ) {}

    quadratic_fitness_parameter( boost::property_tree::ptree & config ) :
        m_scale( 1.0 ) {
        boost::property_tree::ptree fblock;
        fblock = config.get_child( FITNESS_BLOCK_K, fblock );
        
        boost::property_tree::ptree pblock;
        pblock = fblock.get_child( PARAM_K, pblock );

        m_scale = pblock.get< RealType >( SCALE_K, m_scale );
        
        pblock.put( SCALE_K, m_scale );

        fblock.put( NAME_K, QUAD_NAME );
        fblock.put_child( PARAM_K, pblock );

        config.put_child( FITNESS_BLOCK_K, fblock );
    }
};

#endif  // QUADRATIC_FITNESS_PARAMETER_HPP_
