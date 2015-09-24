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
#include "clotho/fitness/random_fitness_generator.hpp"
#include "clotho/fitness/fitness_toolkit.hpp"

const std::string SCALE_K = "scale";
const std::string MU_K = "mu";

random_fitness_generator::random_fitness_generator() :
    m_scale( 1. )
    , m_mu( 1. ) {
    fitness_toolkit::getInstance()->register_tool( this );
}

random_fitness_generator::random_fitness_generator( boost::property_tree::ptree & config ) :
    m_scale(1.)
    , m_mu( 1. ) {
    parseConfig( config );
}

std::shared_ptr< ifitness_generator > random_fitness_generator::create( boost::property_tree::ptree & config ) const {
    std::shared_ptr< ifitness_generator > t( new random_fitness_generator( config ) );
    return t;
}

std::shared_ptr< ifitness > random_fitness_generator::generate( const std::vector< std::vector< double > > & pop_traits ) {
    // theoretical standard deviation:
    // sqrt( 2 * N * mu), where
    //  N - is the haploid sequence count
    //  mu - mutation rate per sequence

    double res = 4.0 * (double)pop_traits.size();
    res *= m_mu;
    res = sqrt( res );  // theoretical standard deviation

    std::shared_ptr< ifitness > t( new result_type( m_scale * res ));
    return t;
}

const std::string random_fitness_generator::name() const {
    return QUAD_NAME;
}

void random_fitness_generator::parseConfig( boost::property_tree::ptree & config ) {
    std::ostringstream oss;
    oss /*<< QUAD_NAME << "."*/ << SCALE_K;
    if( config.get_child_optional( oss.str() ) == boost::none ) {
        config.put( oss.str(), m_scale );
    } else {
        m_scale = config.get< double >( oss.str(), m_scale );

        assert( m_scale > 0.0 );
    }

    oss.str("");
    oss.clear();
    oss /*<< QUAD_NAME << "."*/ << MU_K;
    if( config.get_child_optional( oss.str() ) == boost::none ) {
        config.put( oss.str(), m_mu );
    } else {
        m_mu = config.get< double >( oss.str(), m_mu );
    }
}

void random_fitness_generator::log( std::ostream & out ) const {
    out << "{" << QUAD_NAME << "_gen, " << m_scale << ", " << m_mu << "}\n";
}

random_fitness_generator::~random_fitness_generator() {}

static const random_fitness_generator qfg_reg;
