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
#include "clotho/fitness/quadratic_fitness_generator.hpp"
#include "clotho/fitness/fitness_toolkit.hpp"
#include "clotho/utility/clotho_strings.hpp"
#include <cmath>

quadratic_fitness_generator::quadratic_fitness_generator() :
    quadratic_fitness_parameter< double >()
    , mutation_rate_parameter< double >()
{
    fitness_toolkit::getInstance()->register_tool( this );
}

quadratic_fitness_generator::quadratic_fitness_generator( boost::property_tree::ptree & config ) :
    quadratic_fitness_parameter< double >( config )
    , mutation_rate_parameter< double >( config )
{
//    parseConfig( config );
}

std::shared_ptr< ifitness_generator > quadratic_fitness_generator::create( boost::property_tree::ptree & config ) const {
    std::shared_ptr< ifitness_generator > t( new quadratic_fitness_generator( config ) );
    return t;
}

std::shared_ptr< ifitness > quadratic_fitness_generator::generate( const std::vector< std::vector< real_type > > & pop_traits ) {
    return generate( pop_traits.size() );
}

std::shared_ptr< ifitness > quadratic_fitness_generator::generate( size_t N ) {
    // theoretical standard deviation:
    // sqrt( 2 * N * mu), where
    //  N - is the haploid sequence count
    //  mu - mutation rate per sequence
    //
    double std = 4.0 * (double)N;
    std *= m_mu;
    std = sqrt( std );  // theoretical standard deviation

    std::shared_ptr< ifitness > t( new result_type( m_scale * std ));
    return t;
}

const std::string quadratic_fitness_generator::name() const {
    return QUAD_NAME;
}

//void quadratic_fitness_generator::parseConfig( boost::property_tree::ptree & config ) {
//    boost::property_tree::ptree fblock;
//    fblock = config.get_child( FITNESS_BLOCK_K, fblock );
//    
//    boost::property_tree::ptree pblock;
//    pblock = fblock.get_child( PARAM_K, pblock );
//
//    m_scale = pblock.get< real_type >( SCALE_K, m_scale );
//    
//    pblock.put( SCALE_K, m_scale );
//
//    fblock.put_child( PARAM_K, pblock );
//    config.put_child( FITNESS_BLOCK_K, fblock );
//}

void quadratic_fitness_generator::log( std::ostream & out ) const {
    out << "{" << QUAD_NAME << "_gen, " << m_scale << ", " << m_mu << "}\n";
}

quadratic_fitness_generator::~quadratic_fitness_generator() {}

static const quadratic_fitness_generator qfg_reg;
