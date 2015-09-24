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

const std::string SCALE_K = "scale";
const std::string MU_K = "mu";

quadratic_fitness_generator::quadratic_fitness_generator() :
    m_scale( 1. )
    , m_mu( 1. ) {
    fitness_toolkit::getInstance()->register_tool( this );
}

quadratic_fitness_generator::quadratic_fitness_generator( boost::property_tree::ptree & config ) :
    m_scale(1.)
    , m_mu( 1. ) {
    parseConfig( config );
}

std::shared_ptr< ifitness_generator > quadratic_fitness_generator::create( boost::property_tree::ptree & config ) const {
    std::shared_ptr< ifitness_generator > t( new quadratic_fitness_generator( config ) );
    return t;
}

std::shared_ptr< ifitness > quadratic_fitness_generator::generate( const std::vector< std::vector< double > > & pop_traits ) {
    // theoretical standard deviation:
    // sqrt( 2 * N * mu), where
    //  N - is the haploid sequence count
    //  mu - mutation rate per sequence

    double std = 4.0 * (double)pop_traits.size();
    std *= m_mu;
    std = sqrt( std );  // theoretical standard deviation

    std::shared_ptr< ifitness > t( new result_type( m_scale * std ));
    return t;
}

const std::string quadratic_fitness_generator::name() const {
    return QUAD_NAME;
}

void quadratic_fitness_generator::parseConfig( boost::property_tree::ptree & config ) {
    if( config.get_child_optional( SCALE_K ) == boost::none ) {
        config.put( SCALE_K, m_scale );
    } else {
        m_scale = config.get< double >( SCALE_K, m_scale );

        assert( m_scale > 0.0 );
    }

    if( config.get_child_optional( MU_K ) == boost::none ) {
        config.put( MU_K, m_mu );
    } else {
        m_mu = config.get< double >( MU_K, m_mu );
    }
}

void quadratic_fitness_generator::log( std::ostream & out ) const {
    out << "{" << QUAD_NAME << "_gen, " << m_scale << ", " << m_mu << "}\n";
}

quadratic_fitness_generator::~quadratic_fitness_generator() {}

static const quadratic_fitness_generator qfg_reg;
