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
#include "clotho/fitness/normal_fitness_generator.hpp"
#include "clotho/fitness/fitness_toolkit.hpp"

normal_fitness_generator::normal_fitness_generator() :
    mutation_rate_parameter< double >()
{
    fitness_toolkit::getInstance()->register_tool( this );
}

normal_fitness_generator::normal_fitness_generator( boost::property_tree::ptree & config ) :
    mutation_rate_parameter( config )
{ }

std::shared_ptr< ifitness_generator > normal_fitness_generator::create( boost::property_tree::ptree & config ) const {
    std::shared_ptr< ifitness_generator > t( new normal_fitness_generator( config ) );
    return t;
}

std::shared_ptr< ifitness > normal_fitness_generator::generate( const std::vector< std::vector< double > > & pop_traits ) {
    // theoretical standard deviation:
    // sqrt( 2 * P * N * mu), where
    //  N - is the individual count
    //  P - is the ploidy of an individual = 2.0
    //  mu - mutation rate per sequence
    //
    //  NP = haploid sequence count of the population

    double n = 4.0 * (double)pop_traits.size();

    n *= m_mu;
    n = sqrt(n);    // theoretical standard deviation

    std::shared_ptr< ifitness > r( new result_type( 0., n ) );
    return r;
}

const std::string normal_fitness_generator::name() const {
    return NORM_NAME;
}

void normal_fitness_generator::log( std::ostream & out ) const {
    out << "{" << NORM_NAME << "_gen"
        << ", " << m_mu
        << "}\n";
}

normal_fitness_generator::~normal_fitness_generator() {}

static const normal_fitness_generator nfg_reg;
